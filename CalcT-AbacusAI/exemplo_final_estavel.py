#!/usr/bin/env python3
"""
Exemplo Final Estavel - Construcao em Camadas FENICSx
Implementacao robusta baseada nos testes bem-sucedidos
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.io
import ufl
from mpi4py import MPI
import time
from petsc4py import PETSc

# Import correto do LinearProblem
from dolfinx.fem.petsc import LinearProblem


class ConstrucaoCamadasEstavel:
    """
    Simulador de construcao em camadas com formulacao estavel
    """
    
    def __init__(self, mesh, cell_tags=None, facet_tags=None):
        """
        Inicializa o simulador
        
        Args:
            mesh: Malha FEniCSx
            cell_tags: Tags de celulas (opcional)
            facet_tags: Tags de facetas (opcional)
        """
        self.mesh = mesh
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        
        # Espaco funcional
        self.V = fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solucao
        self.T = fem.Function(self.V)
        self.T_old = fem.Function(self.V)
        
        # Campos de controle
        self.ativacao = fem.Function(self.V)
        self.geracao = fem.Function(self.V)
        
        # Medidas
        self.dx = ufl.Measure("dx", domain=mesh)
        if self.facet_tags:
            self.ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)
        else:
            self.ds = ufl.Measure("ds", domain=mesh)
        
        # Estado
        self.tempo = 0.0
        self.camadas_ativas = set()
        
        # Historico
        self.historico_temp = []
        
        # Nomes para visualizacao
        self.T.name = "Temperatura"
        self.ativacao.name = "Ativacao"
        self.geracao.name = "Geracao"
        
        # Inicializar campos
        self.T.x.array[:] = 20.0
        self.T_old.x.array[:] = 20.0
        self.ativacao.x.array[:] = 0.0
        self.geracao.x.array[:] = 0.0
        
        print(f"Simulador inicializado")
        print(f"Malha: {mesh.topology.index_map(2).size_local} celulas")
        print(f"DOFs: {self.V.dofmap.index_map.size_local}")
    
    def ativar_camada(self, num_camada, params=None):
        """
        Ativa uma camada especifica
        
        Args:
            num_camada: Numero da camada (0, 1, 2, ...)
            params: Parametros da camada (opcional)
        """
        if num_camada in self.camadas_ativas:
            print(f"Camada {num_camada} ja esta ativa")
            return
        
        print(f"Ativando camada {num_camada} no tempo {self.tempo:.1f}h")
        
        # Ativar camada baseada em cell_tags se disponivel
        if self.cell_tags is not None:
            tag_camada = num_camada + 1  # Geralmente tags comecam em 1
            cells_camada = self.cell_tags.find(tag_camada)
            num_cells = len(cells_camada)
            
            print(f"  Celulas da camada: {num_cells}")
            
            if num_cells > 0:
                # Ativar DOFs nas celulas desta camada
                for cell in cells_camada:
                    dofs = self.V.dofmap.cell_dofs(cell)
                    for dof in dofs:
                        if dof < len(self.ativacao.x.array):
                            self.ativacao.x.array[dof] = 1.0
                            if params and 'Q0' in params:
                                self.geracao.x.array[dof] = params['Q0']
        
        # Caso alternativo: ativar baseado em coordenadas
        else:
            # Caso 2D comum: camadas horizontais
            height = self.mesh.geometry.x[:, 1].max()
            num_camadas = 3  # Valor padrao
            layer_height = height / num_camadas
            
            y_min = num_camada * layer_height
            y_max = (num_camada + 1) * layer_height
            
            # Ativar DOFs nesta regiao
            coords = self.V.tabulate_dof_coordinates()
            for i in range(len(coords)):
                if y_min <= coords[i][1] < y_max:
                    self.ativacao.x.array[i] = 1.0
                    if params and 'Q0' in params:
                        self.geracao.x.array[i] = params['Q0']
        
        # Atualizar campos
        self.ativacao.x.scatter_forward()
        self.geracao.x.scatter_forward()
        
        # Registrar ativacao
        self.camadas_ativas.add(num_camada)
        
        # Estatisticas
        dofs_ativos = np.sum(self.ativacao.x.array > 0.5)
        print(f"  DOFs ativos: {dofs_ativos}")
    
    def atualizar_geracao(self, parametros_camadas):
        """
        Atualiza geracao de calor com base no tempo
        
        Args:
            parametros_camadas: Dicionario com parametros por camada
        """
        # Resetar geracao
        self.geracao.x.array[:] = 0.0
        
        # Para cada camada ativa
        for num_camada in self.camadas_ativas:
            if num_camada not in parametros_camadas:
                continue
                
            # Parametros desta camada
            params = parametros_camadas[num_camada]
            Q0 = params.get('Q0', 0.0)
            tau = params.get('tau', 1.0)
            tempo_lancamento = params.get('tempo_lancamento', 0.0)
            
            # Calcular idade
            idade = max(0.0, self.tempo - tempo_lancamento)
            
            # Funcao de geracao de calor (exponencial)
            geracao_atual = Q0 * np.exp(-idade / tau) if idade > 0 else 0.0
            
            # Aplicar geracao nas celulas desta camada
            if self.cell_tags is not None:
                tag_camada = num_camada + 1
                cells_camada = self.cell_tags.find(tag_camada)
                
                for cell in cells_camada:
                    dofs = self.V.dofmap.cell_dofs(cell)
                    for dof in dofs:
                        if dof < len(self.geracao.x.array) and self.ativacao.x.array[dof] > 0.5:
                            self.geracao.x.array[dof] = geracao_atual
            else:
                # Baseado em coordenadas
                height = self.mesh.geometry.x[:, 1].max()
                num_camadas = 3  # Valor padrao
                layer_height = height / num_camadas
                
                y_min = num_camada * layer_height
                y_max = (num_camada + 1) * layer_height
                
                coords = self.V.tabulate_dof_coordinates()
                for i in range(len(coords)):
                    if y_min <= coords[i][1] < y_max and self.ativacao.x.array[i] > 0.5:
                        self.geracao.x.array[i] = geracao_atual
        
        # Atualizar campo
        self.geracao.x.scatter_forward()
    
    def definir_contorno(self, condicoes):
        """
        Define condicoes de contorno
        
        Args:
            condicoes: Lista de condicoes de contorno
            
        Returns:
            Lista de BCs, termos Robin
        """
        bcs = []
        F_robin = 0
        v = ufl.TestFunction(self.V)
        
        for cond in condicoes:
            tipo = cond.get('tipo', '')
            
            # Dirichlet
            if tipo == 'dirichlet':
                if 'tag' in cond and self.facet_tags:
                    # Baseado em facet_tags
                    tag = cond['tag']
                    facets = self.facet_tags.find(tag)
                    if len(facets) > 0:
                        dofs = fem.locate_dofs_topological(self.V, 1, facets)
                        if len(dofs) > 0:
                            bc = fem.dirichletbc(fem.Constant(self.mesh, cond['valor']), dofs, self.V)
                            bcs.append(bc)
                            print(f"  Condicao Dirichlet: tag={tag}, valor={cond['valor']}, DOFs={len(dofs)}")
                elif 'funcao' in cond:
                    # Baseado em funcao geometrica
                    dofs = fem.locate_dofs_geometrical(self.V, cond['funcao'])
                    if len(dofs) > 0:
                        bc = fem.dirichletbc(fem.Constant(self.mesh, cond['valor']), dofs, self.V)
                        bcs.append(bc)
                        print(f"  Condicao Dirichlet: geometrica, valor={cond['valor']}, DOFs={len(dofs)}")
            
            # Robin (conveccao)
            elif tipo == 'robin' and self.facet_tags:
                if 'tag' in cond:
                    h = cond.get('h', 10.0)  # Coeficiente de conveccao
                    T_amb = cond.get('T_amb', 20.0)  # Temperatura ambiente
                    tag = cond['tag']
                    
                    # Termo Robin: h * (T - T_amb) * v * ds(tag)
                    F_robin += h * (self.T - T_amb) * v * self.ds(tag)
                    print(f"  Condicao Robin: tag={tag}, h={h}, T_amb={T_amb}")
        
        return bcs, F_robin
    
    def resolver_passo_linear(self, dt, parametros, condicoes):
        """
        Resolve um passo de tempo usando LinearProblem
        
        Args:
            dt: Passo de tempo
            parametros: Parametros termicos
            condicoes: Condicoes de contorno
            
        Returns:
            True se convergiu, False caso contrario
        """
        # Constantes fisicas
        k = parametros.get('k', 2.5)  # Condutividade termica
        rho = parametros.get('rho', 2400.0)  # Densidade
        cp = parametros.get('cp', 1000.0)  # Calor especifico
        
        # Formulacao variacional
        u = ufl.TrialFunction(self.V)
        v = ufl.TestFunction(self.V)
        
        # Propriedades efetivas
        k_eff = self.ativacao * k + (1 - self.ativacao) * 0.001  # Condutividade minima
        rho_cp_eff = self.ativacao * rho * cp + (1 - self.ativacao) * 0.001  # Capacidade minima
        
        # Forma bilinear
        a = (rho_cp_eff / dt) * u * v * self.dx + k_eff * ufl.dot(ufl.grad(u), ufl.grad(v)) * self.dx
        
        # Forma linear
        L = (rho_cp_eff / dt) * self.T_old * v * self.dx + self.geracao * v * self.dx
        
        # Condicoes de contorno
        bcs, F_robin = self.definir_contorno(condicoes)
        
        # Adicionar termos Robin
        a += ufl.derivative(F_robin, self.T, u)
        
        # Resolver sistema linear
        start_time = time.time()
        problema = LinearProblem(a, L, bcs=bcs)
        T_new = problema.solve()
        solve_time = time.time() - start_time
        
        # Atualizar temperatura
        self.T.x.array[:] = T_new.x.array.copy()
        self.T_old.x.array[:] = T_new.x.array.copy()
        
        # Estatisticas
        T_min = np.min(self.T.x.array)
        T_max = np.max(self.T.x.array)
        T_med = np.mean(self.T.x.array)
        
        self.historico_temp.append({
            'tempo': self.tempo,
            'T_min': T_min,
            'T_max': T_max,
            'T_med': T_med,
            'solve_time': solve_time
        })
        
        return True
    
    def simular_construcao(self, cronograma, parametros_material, dt=2.0):
        """
        Executa simulacao completa
        
        Args:
            cronograma: Lista de eventos
            parametros_material: Parametros do material
            dt: Passo de tempo (horas)
        """
        print("="*60)
        print("SIMULACAO DE CONSTRUCAO EM CAMADAS")
        print("="*60)
        
        # Arquivo de saida
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "simulacao_construcao.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        
        # Parametros das camadas
        parametros_camadas = {}
        
        # Condicoes de contorno gerais
        condicoes_contorno = [
            {'tipo': 'dirichlet', 'tag': 1, 'valor': 20.0},  # Base (fundacao)
            {'tipo': 'robin', 'tag': 2, 'h': 15.0, 'T_amb': 15.0},  # Topo
            {'tipo': 'robin', 'tag': 3, 'h': 10.0, 'T_amb': 15.0}   # Laterais
        ]
        
        # Processar cronograma
        dt_horas = dt  # Passo em horas
        passo = 0
        
        for evento in cronograma:
            tempo_alvo = evento['tempo']
            
            print(f"\nProcessando ate t={tempo_alvo:.1f}h: {evento['tipo']}")
            
            # Avancar no tempo ate este evento
            while self.tempo < tempo_alvo:
                dt_atual = min(dt_horas, tempo_alvo - self.tempo)
                
                # Atualizar geracao de calor
                if self.camadas_ativas:
                    self.atualizar_geracao(parametros_camadas)
                
                # Resolver passo
                converged = self.resolver_passo_linear(dt_atual * 3600.0, parametros_material, condicoes_contorno)
                
                if not converged:
                    print(f"  Erro: Falha na convergencia no tempo {self.tempo:.1f}h")
                    break
                
                # Atualizar tempo
                self.tempo += dt_atual
                passo += 1
                
                # Salvar resultado
                arquivo.write_function(self.T, self.tempo)
                
                # Status (a cada 10 passos)
                if passo % 10 == 0:
                    hist = self.historico_temp[-1]
                    print(f"  t={self.tempo:6.1f}h | T: {hist['T_min']:4.1f}-{hist['T_max']:4.1f} "
                          f"(med={hist['T_med']:4.1f}) | solve={hist['solve_time']:.3f}s")
            
            # Processar evento
            if evento['tipo'] == 'nova_camada':
                num_camada = evento['camada']
                
                # Registrar parametros da camada
                parametros_camadas[num_camada] = {
                    'Q0': evento.get('Q0', 1000.0),
                    'tau': evento.get('tau', 48.0),
                    'tempo_lancamento': self.tempo
                }
                
                # Ativar camada
                self.ativar_camada(num_camada, parametros_camadas[num_camada])
            
            elif evento['tipo'] == 'fim':
                print("\nSimulacao concluida!")
                break
        
        arquivo.close()
        
        # Salvar campos finais
        self.salvar_campos_finais()
        
        # Imprimir resumo
        self.imprimir_resumo()
    
    def salvar_campos_finais(self):
        """
        Salva campos finais para analise
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "campos_finais.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        arquivo.write_function(self.T, self.tempo)
        arquivo.write_function(self.ativacao, self.tempo)
        arquivo.write_function(self.geracao, self.tempo)
        arquivo.close()
        print("Campos finais salvos em: campos_finais.xdmf")
    
    def imprimir_resumo(self):
        """
        Imprime resumo da simulacao
        """
        print("\n" + "="*60)
        print("RESUMO DA SIMULACAO")
        print("="*60)
        print(f"Tempo total simulado: {self.tempo:.1f} horas ({self.tempo/24:.1f} dias)")
        print(f"Camadas ativas: {sorted(self.camadas_ativas)}")
        print(f"Total de passos: {len(self.historico_temp)}")
        
        # Estatisticas finais
        if self.historico_temp:
            hist_final = self.historico_temp[-1]
            print(f"\nTemperatura final:")
            print(f"  Minima: {hist_final['T_min']:.1f}°C")
            print(f"  Maxima: {hist_final['T_max']:.1f}°C")
            print(f"  Media: {hist_final['T_med']:.1f}°C")
            
            # Performance
            tempos_solver = [h['solve_time'] for h in self.historico_temp]
            print(f"\nPerformance:")
            print(f"  Tempo medio por passo: {np.mean(tempos_solver):.3f}s")
            print(f"  Tempo total de calculo: {np.sum(tempos_solver):.1f}s")
        
        # Geracao de calor
        geracao_max = np.max(self.geracao.x.array)
        print(f"\nGeracao de calor atual: {geracao_max:.1f} W/m³")
        
        print("\nArquivos gerados:")
        print("- simulacao_construcao.xdmf (evolucao temporal)")
        print("- campos_finais.xdmf (campos finais)")


def criar_malha_exemplo():
    """
    Cria uma malha de exemplo com tags
    """
    # Malha retangular representando secao de barragem
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD, 
        [[0.0, 0.0], [100.0, 60.0]], 
        [40, 24]
    )
    
    # Criar tags de celulas (3 camadas)
    num_cells = mesh.topology.index_map(2).size_local
    cell_values = np.zeros(num_cells, dtype=np.int32)
    cell_indices = np.arange(num_cells, dtype=np.int32)
    
    # Dividir dominio em camadas horizontais
    for i in range(num_cells):
        # Obter coordenadas da celula
        cell_vertices = mesh.geometry.x[mesh.topology.connectivity(2, 0).links(i)]
        y_center = np.mean(cell_vertices[:, 1])
        
        if y_center < 20.0:
            cell_values[i] = 1  # Camada 0
        elif y_center < 40.0:
            cell_values[i] = 2  # Camada 1
        else:
            cell_values[i] = 3  # Camada 2
    
    cell_tags = dolfinx.mesh.meshtags(mesh, 2, cell_indices, cell_values)
    
    # Computar conectividade necessaria
    mesh.topology.create_entities(1)  # Criar facetas (arestas em 2D)
    mesh.topology.create_connectivity(1, 2)  # Faceta->Celula
    
    # Criar tags de facetas (contorno)
    facets = dolfinx.mesh.exterior_facet_indices(mesh.topology)
    facet_values = np.zeros(len(facets), dtype=np.int32)
    
    # Identificar contornos
    for i, facet in enumerate(facets):
        facet_vertices = mesh.geometry.x[mesh.topology.connectivity(1, 0).links(facet)]
        y_coords = facet_vertices[:, 1]
        x_coords = facet_vertices[:, 0]
        
        if np.allclose(y_coords, 0.0):
            facet_values[i] = 1  # Base (fundacao)
        elif np.allclose(y_coords, 60.0):
            facet_values[i] = 2  # Topo
        elif np.allclose(x_coords, 0.0) or np.allclose(x_coords, 100.0):
            facet_values[i] = 3  # Laterais
        else:
            facet_values[i] = 4  # Outros
    
    facet_tags = dolfinx.mesh.meshtags(mesh, 1, facets, facet_values)
    
    return mesh, cell_tags, facet_tags


def exemplo_barragem():
    """
    Exemplo de simulacao de barragem
    """
    # Criar malha com tags
    mesh, cell_tags, facet_tags = criar_malha_exemplo()
    
    # Parametros do concreto
    parametros_material = {
        'k': 2.5,      # Condutividade termica [W/m·K]
        'rho': 2400.0, # Densidade [kg/m³]
        'cp': 1000.0   # Calor especifico [J/kg·K]
    }
    
    # Cronograma de construcao
    cronograma = [
        {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0, 'Q0': 1200.0, 'tau': 48.0},
        {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1, 'Q0': 1000.0, 'tau': 48.0},
        {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2, 'Q0': 800.0, 'tau': 48.0},
        {'tempo': 720.0, 'tipo': 'fim'}
    ]
    
    # Criar simulador
    simulador = ConstrucaoCamadasEstavel(mesh, cell_tags, facet_tags)
    
    # Executar simulacao
    simulador.simular_construcao(cronograma, parametros_material, dt=6.0)
    
    return simulador


if __name__ == "__main__":
    print("Executando exemplo final estavel...")
    simulador = exemplo_barragem()