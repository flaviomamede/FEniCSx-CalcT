#!/usr/bin/env python3
"""
Exemplo Completo - Construção em Camadas FENICSx
Demonstra o uso integrado de todas as funcionalidades
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.io
import dolfinx.nls.petsc
import ufl
from mpi4py import MPI
import time


def criar_malha_exemplo():
    """
    Cria malha de exemplo para demonstração
    """
    # Malha retangular representando seção de barragem
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD, 
        [[0.0, 0.0], [100.0, 60.0]], 
        [40, 24]
    )
    
    # Criar tags de células (3 camadas)
    num_cells = mesh.topology.index_map(2).size_local
    cell_values = np.zeros(num_cells, dtype=np.int32)
    
    # Dividir domínio em camadas horizontais
    for i in range(num_cells):
        # Obter coordenadas da célula
        cell_vertices = mesh.geometry.x[mesh.topology.connectivity(2, 0).links(i)]
        y_center = np.mean(cell_vertices[:, 1])
        
        if y_center < 20.0:
            cell_values[i] = 1  # Camada 0
        elif y_center < 40.0:
            cell_values[i] = 2  # Camada 1
        else:
            cell_values[i] = 3  # Camada 2
    
    cell_tags = dolfinx.mesh.meshtags(mesh, 2, np.arange(num_cells), cell_values)
    
    # Criar tags de facetas (contorno)
    mesh.topology.create_connectivity(1, 2)
    boundary_facets = dolfinx.mesh.exterior_facet_indices(mesh.topology)
    facet_values = np.zeros(len(boundary_facets), dtype=np.int32)
    
    # Identificar contornos
    for i, facet in enumerate(boundary_facets):
        facet_vertices = mesh.geometry.x[mesh.topology.connectivity(1, 0).links(facet)]
        y_coords = facet_vertices[:, 1]
        x_coords = facet_vertices[:, 0]
        
        if np.allclose(y_coords, 0.0):
            facet_values[i] = 1  # Base (fundação)
        elif np.allclose(y_coords, 60.0):
            facet_values[i] = 2  # Topo
        elif np.allclose(x_coords, 0.0) or np.allclose(x_coords, 100.0):
            facet_values[i] = 3  # Laterais
        else:
            facet_values[i] = 4  # Outros
    
    facet_tags = dolfinx.mesh.meshtags(mesh, 1, boundary_facets, facet_values)
    
    return mesh, cell_tags, facet_tags


class BarragemConstrucaoCamadas:
    """
    Simulador específico para barragem com construção em camadas
    """
    
    def __init__(self, mesh, cell_tags, facet_tags):
        self.mesh = mesh
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        
        # Espaço funcional
        self.V = fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solução
        self.T = fem.Function(self.V)
        self.T_old = fem.Function(self.V)
        
        # Campos de controle
        self.ativacao = fem.Function(self.V)
        self.geracao = fem.Function(self.V)
        
        # Medidas
        self.dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags)
        self.ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)
        
        # Estado
        self.tempo = 0.0
        self.camadas_ativas = set()
        self.solver = None
        
        # Histórico de temperaturas
        self.historico_temp = []
        
        # Inicializar
        self.T.x.array[:] = 20.0  # Temperatura inicial
        self.T_old.x.array[:] = 20.0
        self.ativacao.x.array[:] = 0.0
        self.geracao.x.array[:] = 0.0
        
        print(f"Simulador de barragem inicializado")
        print(f"Malha: {mesh.topology.index_map(2).size_local} células")
        print(f"DOFs: {self.V.dofmap.index_map.size_local}")
    
    def ativar_camada(self, num_camada, parametros_camada):
        """
        Ativa uma camada com parâmetros específicos
        """
        if num_camada in self.camadas_ativas:
            return
        
        print(f"Ativando camada {num_camada} no tempo {self.tempo:.1f}h")
        
        # Adicionar à lista
        self.camadas_ativas.add(num_camada)
        
        # Tag da camada (1, 2, 3, ...)
        tag_camada = num_camada + 1
        
        # Encontrar células da camada
        cells_camada = self.cell_tags.find(tag_camada)
        
        print(f"  Células da camada: {len(cells_camada)}")
        
        # Ativar células
        for cell in cells_camada:
            dofs = self.V.dofmap.cell_dofs(cell)
            self.ativacao.x.array[dofs] = 1.0
        
        # Atualizar campos
        self.ativacao.x.scatter_forward()
        self.atualizar_geracao(parametros_camada)
        
        # Invalidar solver
        self.solver = None
        
        # Estatísticas
        dofs_ativos = np.sum(self.ativacao.x.array > 0.5)
        print(f"  DOFs ativos: {dofs_ativos}")
    
    def atualizar_geracao(self, parametros_camada):
        """
        Atualiza geração de calor para todas as camadas ativas
        """
        # Resetar geração
        self.geracao.x.array[:] = 0.0
        
        # Processar cada camada ativa
        for num_camada in self.camadas_ativas:
            tag_camada = num_camada + 1
            cells_camada = self.cell_tags.find(tag_camada)
            
            # Parâmetros da camada
            tempo_lancamento = parametros_camada[num_camada]['tempo_lancamento']
            Q0 = parametros_camada[num_camada]['Q0']
            tau = parametros_camada[num_camada]['tau']
            
            # Calcular idade
            idade = max(0.0, self.tempo - tempo_lancamento)
            
            if idade <= 0:
                continue
            
            # Calcular geração
            geracao_atual = Q0 * np.exp(-idade / tau)
            
            # Aplicar nas células da camada
            for cell in cells_camada:
                dofs = self.V.dofmap.cell_dofs(cell)
                self.geracao.x.array[dofs] = geracao_atual
        
        self.geracao.x.scatter_forward()
    
    def criar_formulacao(self, dt, parametros_material):
        """
        Cria formulação variacional
        """
        v = ufl.TestFunction(self.V)
        
        # Propriedades efetivas
        k = self.ativacao * parametros_material['k']
        rho_cp = self.ativacao * parametros_material['rho'] * parametros_material['cp']
        
        # Formulação variacional
        F = (rho_cp * (self.T - self.T_old) / dt) * v * self.dx
        F += k * ufl.dot(ufl.grad(self.T), ufl.grad(v)) * self.dx
        F -= self.geracao * v * self.dx
        
        return F
    
    def definir_contorno(self, condicoes_contorno):
        """
        Define condições de contorno
        """
        bcs = []
        F_robin = 0
        
        for cond in condicoes_contorno:
            if cond['tipo'] == 'dirichlet':
                facets = self.facet_tags.find(cond['tag'])
                if len(facets) > 0:
                    dofs = fem.locate_dofs_topological(self.V, 1, facets)
                    if len(dofs) > 0:
                        bc = fem.dirichletbc(fem.Constant(self.mesh, cond['valor']), dofs, self.V)
                        bcs.append(bc)
                        print(f"  Condição Dirichlet: tag={cond['tag']}, valor={cond['valor']}, DOFs={len(dofs)}")
            
            elif cond['tipo'] == 'robin':
                v = ufl.TestFunction(self.V)
                h = cond['h']
                T_amb = cond['T_amb']
                F_robin += h * (self.T - T_amb) * v * self.ds(cond['tag'])
                print(f"  Condição Robin: tag={cond['tag']}, h={h}, T_amb={T_amb}")
        
        return bcs, F_robin
    
    def resolver_passo(self, dt, parametros_material, condicoes_contorno):
        """
        Resolve um passo de tempo
        """
        # Criar formulação
        F = self.criar_formulacao(dt, parametros_material)
        
        # Condições de contorno
        bcs, F_robin = self.definir_contorno(condicoes_contorno)
        F += F_robin
        
        # Criar/atualizar problema
        problem = fem.petsc.NonlinearProblem(F, self.T, bcs)
        
        # Configurar solver
        if self.solver is None:
            self.solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, problem)
            self.solver.rtol = 1e-6
            self.solver.atol = 1e-10
            self.solver.max_it = 25
            print(f"  Solver reconfigurado com {len(bcs)} condições Dirichlet")
        
        # Resolver
        start_time = time.time()
        n_iter, converged = self.solver.solve(self.T)
        solve_time = time.time() - start_time
        
        if not converged:
            print(f"  Aviso: Não convergiu em {n_iter} iterações")
        
        # Atualizar solução anterior
        self.T_old.x.array[:] = self.T.x.array.copy()
        
        # Estatísticas
        T_min = np.min(self.T.x.array)
        T_max = np.max(self.T.x.array)
        T_med = np.mean(self.T.x.array)
        
        self.historico_temp.append({
            'tempo': self.tempo,
            'T_min': T_min,
            'T_max': T_max,
            'T_med': T_med,
            'n_iter': n_iter,
            'solve_time': solve_time
        })
        
        return converged
    
    def simular_barragem(self, cronograma, parametros_material, dt=2.0):
        """
        Executa simulação completa da barragem
        """
        print("="*60)
        print("SIMULAÇÃO DE BARRAGEM - CONSTRUÇÃO EM CAMADAS")
        print("="*60)
        
        # Arquivo de saída
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "barragem_simulacao.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        
        # Configurar nomes
        self.T.name = "Temperatura"
        self.ativacao.name = "Ativacao"
        self.geracao.name = "Geracao"
        
        # Condições de contorno
        condicoes = [
            {'tipo': 'dirichlet', 'tag': 1, 'valor': 20.0},  # Base (fundação)
            {'tipo': 'robin', 'tag': 2, 'h': 15.0, 'T_amb': 15.0},  # Topo
            {'tipo': 'robin', 'tag': 3, 'h': 10.0, 'T_amb': 15.0},  # Laterais
        ]
        
        # Parâmetros das camadas
        parametros_camadas = {
            0: {'tempo_lancamento': 0.0, 'Q0': 1200.0, 'tau': 48.0},
            1: {'tempo_lancamento': 168.0, 'Q0': 1000.0, 'tau': 48.0},
            2: {'tempo_lancamento': 336.0, 'Q0': 800.0, 'tau': 48.0},
        }
        
        # Processar cronograma
        for evento in cronograma:
            tempo_alvo = evento['tempo']
            
            print(f"\nProcessando até t={tempo_alvo:.1f}h: {evento['tipo']}")
            
            # Avançar no tempo
            passo = 0
            while self.tempo < tempo_alvo:
                dt_atual = min(dt, tempo_alvo - self.tempo)
                
                # Atualizar geração (varia com o tempo)
                if self.camadas_ativas:
                    self.atualizar_geracao(parametros_camadas)
                
                # Resolver passo
                converged = self.resolver_passo(dt_atual, parametros_material, condicoes)
                
                if not converged:
                    print(f"  Erro: Falha na convergência no tempo {self.tempo:.1f}h")
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
                          f"(med={hist['T_med']:4.1f}) | iter={hist['n_iter']:2d} | "
                          f"solve={hist['solve_time']:5.3f}s")
            
            # Processar evento
            if evento['tipo'] == 'nova_camada':
                self.ativar_camada(evento['camada'], parametros_camadas)
            elif evento['tipo'] == 'fim':
                print("\nSimulação concluída!")
                break
        
        arquivo.close()
        
        # Salvar campos finais
        self.salvar_campos_finais()
        
        # Resumo final
        self.imprimir_resumo()
    
    def salvar_campos_finais(self):
        """
        Salva campos finais para análise
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "barragem_campos_finais.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        arquivo.write_function(self.T, self.tempo)
        arquivo.write_function(self.ativacao, self.tempo)
        arquivo.write_function(self.geracao, self.tempo)
        arquivo.close()
        print("Campos finais salvos em: barragem_campos_finais.xdmf")
    
    def imprimir_resumo(self):
        """
        Imprime resumo da simulação
        """
        print("\n" + "="*60)
        print("RESUMO DA SIMULAÇÃO")
        print("="*60)
        print(f"Tempo total simulado: {self.tempo:.1f} horas ({self.tempo/24:.1f} dias)")
        print(f"Camadas ativas: {sorted(self.camadas_ativas)}")
        print(f"Total de passos: {len(self.historico_temp)}")
        
        # Estatísticas finais
        hist_final = self.historico_temp[-1]
        print(f"\nTemperatura final:")
        print(f"  Mínima: {hist_final['T_min']:.1f}°C")
        print(f"  Máxima: {hist_final['T_max']:.1f}°C")
        print(f"  Média: {hist_final['T_med']:.1f}°C")
        
        # Performance
        tempos_solver = [h['solve_time'] for h in self.historico_temp]
        print(f"\nPerformance:")
        print(f"  Tempo médio por passo: {np.mean(tempos_solver):.3f}s")
        print(f"  Tempo total de cálculo: {np.sum(tempos_solver):.1f}s")
        
        # Geração de calor
        geracao_max = np.max(self.geracao.x.array)
        print(f"\nGeração de calor atual: {geracao_max:.1f} W/m³")
        
        print("\nArquivos gerados:")
        print("- barragem_simulacao.xdmf (evolução temporal)")
        print("- barragem_campos_finais.xdmf (campos finais)")


def main():
    """
    Função principal do exemplo
    """
    # Criar malha
    mesh, cell_tags, facet_tags = criar_malha_exemplo()
    
    # Parâmetros do concreto
    parametros_material = {
        'k': 2.5,      # Condutividade térmica [W/m·K]
        'rho': 2400.0, # Densidade [kg/m³]
        'cp': 1000.0   # Calor específico [J/kg·K]
    }
    
    # Cronograma de construção
    cronograma = [
        {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},      # Fundação
        {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},    # 1ª elevação (7 dias)
        {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2},    # 2ª elevação (14 dias)
        {'tempo': 720.0, 'tipo': 'fim'}                          # Fim (30 dias)
    ]
    
    # Criar simulador
    barragem = BarragemConstrucaoCamadas(mesh, cell_tags, facet_tags)
    
    # Executar simulação
    barragem.simular_barragem(cronograma, parametros_material, dt=2.0)


if __name__ == "__main__":
    main()