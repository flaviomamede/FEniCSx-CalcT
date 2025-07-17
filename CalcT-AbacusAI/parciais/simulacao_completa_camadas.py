"""
Exemplo completo de simulação de construção em camadas
Integração completa com solver e loop temporal
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.fem.petsc
import dolfinx.io
import dolfinx.nls.petsc
import ufl
from mpi4py import MPI
from petsc4py import PETSc
import time


class SimuladorConstrucaoCamadas:
    """
    Simulador completo para construção em camadas
    """
    
    def __init__(self, mesh, cell_tags, facet_tags, camadas_info, parametros_termicos):
        self.mesh = mesh
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        self.camadas_info = camadas_info
        self.parametros_termicos = parametros_termicos
        
        # Configurar espaço funcional
        self.V = fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solução
        self.T = fem.Function(self.V)
        self.T_old = fem.Function(self.V)
        
        # Campos de ativação
        self.ativacao_total = fem.Function(self.V)
        self.geracao_total = fem.Function(self.V)
        
        # Medidas de integração
        self.dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags)
        self.ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)
        
        # Estado da simulação
        self.tempo_atual = 0.0
        self.camadas_ativas = set()
        
        # Solver
        self.solver = None
        self.problem = None
        
        # Inicializar temperatura
        self.T.x.array[:] = 20.0
        self.T_old.x.array[:] = 20.0
        
        # Inicializar campos
        self.ativacao_total.x.array[:] = 0.0
        self.geracao_total.x.array[:] = 0.0
        
        print(f"Simulador inicializado com {len(camadas_info)} camadas")
    
    def ativar_camada(self, num_camada):
        """
        Ativa uma nova camada
        """
        if num_camada in self.camadas_ativas:
            return
        
        print(f"Ativando camada {num_camada} no tempo {self.tempo_atual:.1f}")
        
        # Adicionar à lista de camadas ativas
        self.camadas_ativas.add(num_camada)
        
        # Atualizar campos
        self.atualizar_campos_ativacao()
        
        # Invalidar solver (precisa reconfigurar)
        self.solver = None
        self.problem = None
    
    def atualizar_campos_ativacao(self):
        """
        Atualiza campos de ativação e geração
        """
        # Resetar campos
        self.ativacao_total.x.array[:] = 0.0
        self.geracao_total.x.array[:] = 0.0
        
        # Processar cada camada ativa
        for num_camada in self.camadas_ativas:
            camada_info = self.camadas_info[num_camada]
            
            # Obter células da camada
            cells_camada = self.obter_cells_camada(num_camada)
            
            if len(cells_camada) == 0:
                continue
            
            # Calcular idade da camada
            idade = max(0.0, self.tempo_atual - camada_info['tempo_lancamento'])
            
            # Calcular geração de calor
            Q0 = camada_info['material']['Q0']
            tau = camada_info['material']['tau']
            geracao = Q0 * np.exp(-idade / tau) if idade > 0 else 0.0
            
            # Aplicar nas células da camada
            for cell in cells_camada:
                dofs_cell = self.V.dofmap.cell_dofs(cell)
                self.ativacao_total.x.array[dofs_cell] = 1.0
                self.geracao_total.x.array[dofs_cell] = geracao
        
        # Atualizar
        self.ativacao_total.x.scatter_forward()
        self.geracao_total.x.scatter_forward()
    
    def obter_cells_camada(self, num_camada):
        """
        Retorna células da camada
        """
        tag_camada = num_camada + 1
        cells = []
        
        if hasattr(self.cell_tags, 'find'):
            # Usar método find se disponível
            cells = self.cell_tags.find(tag_camada)
        else:
            # Busca manual
            for i, tag in enumerate(self.cell_tags.values):
                if tag == tag_camada:
                    cells.append(i)
        
        return np.array(cells)
    
    def criar_formulacao_variacional(self, dt):
        """
        Cria formulação variacional
        """
        # Função teste
        v = ufl.TestFunction(self.V)
        
        # Propriedades efetivas
        k_base = self.parametros_termicos['condutividade']
        rho_base = self.parametros_termicos['densidade']
        cp_base = self.parametros_termicos['calor_especifico']
        
        # Aplicar ativação
        k_eff = self.ativacao_total * k_base
        rho_cp_eff = self.ativacao_total * rho_base * cp_base
        
        # Formulação variacional
        # Termo transiente
        F = (rho_cp_eff * (self.T - self.T_old) / dt) * v * self.dx
        
        # Termo de difusão
        F += k_eff * ufl.dot(ufl.grad(self.T), ufl.grad(v)) * self.dx
        
        # Termo de geração
        F -= self.geracao_total * v * self.dx
        
        # Condições de contorno de Robin (convecção)
        for num_camada in self.camadas_ativas:
            camada_info = self.camadas_info[num_camada]
            
            if 'robin' in camada_info:
                for bc_info in camada_info['robin']:
                    tag = bc_info['tag']
                    h = bc_info['h']
                    T_amb = bc_info['T_amb']
                    
                    # Aplicar condição de Robin
                    F += h * (self.T - T_amb) * v * self.ds(tag)
        
        return F
    
    def definir_condicoes_dirichlet(self):
        """
        Define condições de contorno de Dirichlet
        """
        bcs = []
        
        for num_camada in self.camadas_ativas:
            camada_info = self.camadas_info[num_camada]
            
            if 'dirichlet' in camada_info:
                for bc_info in camada_info['dirichlet']:
                    tag = bc_info['tag']
                    valor = bc_info['valor']
                    
                    # Encontrar facetas com a tag
                    facets = self.facet_tags.find(tag)
                    
                    if len(facets) > 0:
                        # Localizar DOFs
                        dofs = fem.locate_dofs_topological(
                            self.V, self.mesh.topology.dim - 1, facets
                        )
                        
                        # Criar condição de contorno
                        bc_value = fem.Constant(self.mesh, valor)
                        bc = fem.dirichletbc(bc_value, dofs, self.V)
                        bcs.append(bc)
        
        return bcs
    
    def configurar_solver(self, dt):
        """
        Configura solver não-linear
        """
        # Criar formulação
        F = self.criar_formulacao_variacional(dt)
        
        # Condições de contorno
        bcs = self.definir_condicoes_dirichlet()
        
        # Criar problema não-linear
        self.problem = fem.petsc.NonlinearProblem(F, self.T, bcs)
        
        # Configurar solver
        self.solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, self.problem)
        
        # Configurações de convergência
        self.solver.convergence_criterion = "incremental"
        self.solver.rtol = 1e-6
        self.solver.atol = 1e-10
        self.solver.max_it = 25
        
        # Configurar solver linear
        ksp = self.solver.krylov_solver
        opts = PETSc.Options()
        opts["ksp_type"] = "preonly"
        opts["pc_type"] = "lu"
        opts["pc_factor_mat_solver_type"] = "mumps"
        ksp.setFromOptions()
        
        print(f"Solver configurado com {len(bcs)} condições de Dirichlet")
    
    def resolver_passo_tempo(self, dt):
        """
        Resolve um passo de tempo
        """
        # Atualizar campos (geração varia com o tempo)
        self.atualizar_campos_ativacao()
        
        # Configurar solver se necessário
        if self.solver is None:
            self.configurar_solver(dt)
        
        # Resolver sistema
        n_iter, converged = self.solver.solve(self.T)
        
        if not converged:
            print(f"Aviso: Convergência não alcançada em {n_iter} iterações")
            return False
        
        # Atualizar solução anterior
        self.T_old.x.array[:] = self.T.x.array.copy()
        
        return True
    
    def simular(self, cronograma, dt_base=1.0, arquivo_saida="simulacao.xdmf"):
        """
        Executa simulação completa
        """
        print("Iniciando simulação...")
        print(f"Cronograma: {len(cronograma)} eventos")
        
        # Configurar saída
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, arquivo_saida, "w")
        arquivo.write_mesh(self.mesh)
        
        # Configurar temperatura para saída
        self.T.name = "Temperatura"
        
        # Processar cronograma
        for evento in cronograma:
            tempo_alvo = evento['tempo']
            
            print(f"Processando evento: {evento['tipo']} no tempo {tempo_alvo}")
            
            # Avançar no tempo até o evento
            while self.tempo_atual < tempo_alvo:
                dt = min(dt_base, tempo_alvo - self.tempo_atual)
                
                # Resolver passo de tempo
                start_time = time.time()
                sucesso = self.resolver_passo_tempo(dt)
                solve_time = time.time() - start_time
                
                if not sucesso:
                    print(f"Erro: Falha na convergência no tempo {self.tempo_atual}")
                    break
                
                # Atualizar tempo
                self.tempo_atual += dt
                
                # Salvar resultado
                arquivo.write_function(self.T, self.tempo_atual)
                
                # Estatísticas
                T_min = np.min(self.T.x.array)
                T_max = np.max(self.T.x.array)
                T_med = np.mean(self.T.x.array)
                geracao_max = np.max(self.geracao_total.x.array)
                
                print(f"t={self.tempo_atual:6.1f}h | "
                      f"T: {T_min:5.1f}-{T_max:5.1f} (med={T_med:5.1f}) | "
                      f"Q_max={geracao_max:6.1f} | "
                      f"solve={solve_time:5.3f}s")
            
            # Processar evento
            if evento['tipo'] == 'nova_camada':
                self.ativar_camada(evento['camada'])
            elif evento['tipo'] == 'fim_simulacao':
                print("Simulação concluída")
                break
        
        arquivo.close()
        
        # Salvar campos finais
        self.salvar_campos_debug("campos_finais.xdmf")
        
        print(f"Simulação finalizada. Tempo total: {self.tempo_atual:.1f}h")
    
    def salvar_campos_debug(self, nome_arquivo):
        """
        Salva campos para debug
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, nome_arquivo, "w")
        arquivo.write_mesh(self.mesh)
        
        # Configurar nomes
        self.T.name = "Temperatura"
        self.ativacao_total.name = "Ativacao"
        self.geracao_total.name = "Geracao"
        
        # Salvar campos
        arquivo.write_function(self.T, self.tempo_atual)
        arquivo.write_function(self.ativacao_total, self.tempo_atual)
        arquivo.write_function(self.geracao_total, self.tempo_atual)
        
        arquivo.close()
        print(f"Campos salvos em {nome_arquivo}")


def exemplo_simulacao_completa():
    """
    Exemplo completo de simulação
    """
    # Criar malha simples para teste
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [100.0, 60.0]],
        [50, 30]
    )
    
    # Criar tags manuais (simulando camadas)
    num_cells = mesh.topology.index_map(mesh.topology.dim).size_local
    cell_values = np.zeros(num_cells, dtype=np.int32)
    
    # Dividir em 3 camadas (bottom to top)
    for i in range(num_cells):
        cell = mesh.geometry.x[mesh.topology.connectivity(2, 0).links(i)]
        y_center = np.mean(cell[:, 1])
        
        if y_center < 20.0:
            cell_values[i] = 1  # Camada 0
        elif y_center < 40.0:
            cell_values[i] = 2  # Camada 1
        else:
            cell_values[i] = 3  # Camada 2
    
    cell_tags = dolfinx.mesh.meshtags(mesh, mesh.topology.dim, 
                                      np.arange(num_cells), cell_values)
    
    # Criar tags de facetas (contorno)
    mesh.topology.create_connectivity(mesh.topology.dim - 1, mesh.topology.dim)
    num_facets = mesh.topology.index_map(mesh.topology.dim - 1).size_local
    facet_values = np.zeros(num_facets, dtype=np.int32)
    
    # Identificar facetas do contorno
    boundary_facets = dolfinx.mesh.exterior_facet_indices(mesh.topology)
    
    for facet in boundary_facets:
        facet_coords = mesh.geometry.x[mesh.topology.connectivity(1, 0).links(facet)]
        y_coords = facet_coords[:, 1]
        
        if np.allclose(y_coords, 0.0):
            facet_values[facet] = 1  # Base
        elif np.allclose(y_coords, 60.0):
            facet_values[facet] = 2  # Topo
        else:
            facet_values[facet] = 3  # Laterais
    
    facet_tags = dolfinx.mesh.meshtags(mesh, mesh.topology.dim - 1,
                                       np.arange(num_facets), facet_values)
    
    # Parâmetros da simulação
    parametros_termicos = {
        'condutividade': 2.5,
        'densidade': 2400.0,
        'calor_especifico': 1000.0
    }
    
    # Configuração das camadas
    camadas_info = [
        {
            'tempo_lancamento': 0.0,
            'material': {'Q0': 1000.0, 'tau': 48.0},
            'dirichlet': [{'tag': 1, 'valor': 20.0}],  # Base fixa
            'robin': [{'tag': 3, 'h': 10.0, 'T_amb': 15.0}]  # Laterais
        },
        {
            'tempo_lancamento': 168.0,
            'material': {'Q0': 1000.0, 'tau': 48.0},
            'robin': [{'tag': 3, 'h': 10.0, 'T_amb': 15.0}]
        },
        {
            'tempo_lancamento': 336.0,
            'material': {'Q0': 1000.0, 'tau': 48.0},
            'robin': [
                {'tag': 3, 'h': 10.0, 'T_amb': 15.0},  # Laterais
                {'tag': 2, 'h': 15.0, 'T_amb': 15.0}   # Topo
            ]
        }
    ]
    
    # Cronograma
    cronograma = [
        {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},
        {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},
        {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2},
        {'tempo': 720.0, 'tipo': 'fim_simulacao'}
    ]
    
    # Criar simulador
    simulador = SimuladorConstrucaoCamadas(
        mesh, cell_tags, facet_tags, camadas_info, parametros_termicos
    )
    
    # Executar simulação
    simulador.simular(cronograma, dt_base=2.0, arquivo_saida="construcao_camadas.xdmf")
    
    return simulador


if __name__ == "__main__":
    print("="*60)
    print("SIMULAÇÃO DE CONSTRUÇÃO EM CAMADAS - FENICSx")
    print("="*60)
    
    simulador = exemplo_simulacao_completa()
    
    print("\nSimulação concluída!")
    print("Arquivos gerados:")
    print("- construcao_camadas.xdmf (resultados temporais)")
    print("- campos_finais.xdmf (campos finais)")