"""
Construção em Camadas com Ativação/Desativação de Elementos - FENICSx 0.9.0
Análise Térmica com Geração Interna de Calor Exotérmica

Estratégia:
1. Malha global única com tags por camada
2. Submeshes ativas controladas por cell_tags
3. Campos de ativação para propriedades
4. Condições de contorno dinâmicas
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.io
import dolfinx.mesh
import dolfinx.plot
import ufl
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx.nls.petsc
import dolfinx.fem.petsc


class ConstrucaoCamadas:
    """
    Classe para gerenciar construção em camadas com ativação/desativação
    """
    
    def __init__(self, mesh, camadas_info, parametros_termicos):
        self.mesh = mesh
        self.camadas_info = camadas_info  # Dict com info de cada camada
        self.parametros_termicos = parametros_termicos
        self.dim = mesh.topology.dim
        self.camada_atual = 0
        self.tempo_atual = 0.0
        
        # Inicializar espaços funcionais
        self.V = fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solução
        self.T = fem.Function(self.V)  # Temperatura atual
        self.T_old = fem.Function(self.V)  # Temperatura anterior
        
        # Campos de ativação
        self.campo_ativacao = fem.Function(self.V)
        self.campo_geracao = fem.Function(self.V)
        
        # Tags de células e facetas
        self.cell_tags = dolfinx.mesh.meshtags_from_entities(
            mesh, mesh.topology.dim, 
            np.arange(mesh.topology.index_map(mesh.topology.dim).size_local),
            np.zeros(mesh.topology.index_map(mesh.topology.dim).size_local, dtype=np.int32)
        )
        
        self.facet_tags = dolfinx.mesh.meshtags_from_entities(
            mesh, mesh.topology.dim - 1,
            np.arange(mesh.topology.index_map(mesh.topology.dim - 1).size_local),
            np.zeros(mesh.topology.index_map(mesh.topology.dim - 1).size_local, dtype=np.int32)
        )
        
        # Medidas de integração
        self.dx = ufl.Measure("dx", domain=mesh, subdomain_data=self.cell_tags)
        self.ds = ufl.Measure("ds", domain=mesh, subdomain_data=self.facet_tags)
        
        # Solver setup
        self.solver = None
        self.problem = None
        
    def atualizar_ativacao_camada(self, nova_camada):
        """
        Atualiza os campos de ativação para nova camada
        """
        print(f"Ativando camada {nova_camada}")
        
        # Atualizar campo de ativação
        campo_ativacao_array = self.campo_ativacao.x.array
        campo_geracao_array = self.campo_geracao.x.array
        
        for i, camada in enumerate(self.camadas_info):
            if i <= nova_camada:
                # Camada ativa
                cells_camada = self.obter_cells_camada(i)
                
                # Ativar elementos da camada
                for cell in cells_camada:
                    dofs_cell = self.V.dofmap.cell_dofs(cell)
                    campo_ativacao_array[dofs_cell] = 1.0
                    
                    # Definir geração de calor baseada na idade da camada
                    idade_camada = self.tempo_atual - camada['tempo_lancamento']
                    geracao = self.calcular_geracao_calor(idade_camada, camada['material'])
                    campo_geracao_array[dofs_cell] = geracao
            else:
                # Camada inativa
                cells_camada = self.obter_cells_camada(i)
                for cell in cells_camada:
                    dofs_cell = self.V.dofmap.cell_dofs(cell)
                    campo_ativacao_array[dofs_cell] = 0.0
                    campo_geracao_array[dofs_cell] = 0.0
        
        # Atualizar os campos
        self.campo_ativacao.x.scatter_forward()
        self.campo_geracao.x.scatter_forward()
        
        self.camada_atual = nova_camada
        
    def obter_cells_camada(self, num_camada):
        """
        Retorna as células que pertencem à camada especificada
        """
        # Implementar lógica baseada na geometria e tags
        # Exemplo: assumindo que as tags foram definidas na criação da malha
        cells = []
        for i, tag in enumerate(self.cell_tags.values):
            if tag == num_camada:
                cells.append(i)
        return cells
    
    def calcular_geracao_calor(self, idade, material):
        """
        Calcula geração de calor exotérmica baseada na idade da camada
        """
        # Exemplo: lei exponencial para geração de calor
        Q0 = material.get('Q0', 1000.0)  # Geração inicial
        tau = material.get('tau', 24.0)  # Constante de tempo (horas)
        
        if idade <= 0:
            return 0.0
            
        return Q0 * np.exp(-idade / tau)
    
    def criar_submesh_ativa(self):
        """
        Cria submesh apenas com células ativas
        """
        # Obter células ativas
        cells_ativas = []
        for i in range(self.camada_atual + 1):
            cells_ativas.extend(self.obter_cells_camada(i))
        
        if not cells_ativas:
            return None
            
        # Criar submesh
        submesh, cell_map, vertex_map, geom_map = dolfinx.mesh.create_submesh(
            self.mesh, self.mesh.topology.dim, cells_ativas
        )
        
        return submesh, cell_map, vertex_map, geom_map
    
    def definir_propriedades_material(self):
        """
        Define propriedades do material com ativação
        """
        # Condutividade térmica efetiva
        k_ref = self.parametros_termicos['condutividade']
        k_eff = self.campo_ativacao * k_ref
        
        # Capacidade térmica efetiva
        rho_cp_ref = self.parametros_termicos['densidade'] * self.parametros_termicos['calor_especifico']
        rho_cp_eff = self.campo_ativacao * rho_cp_ref
        
        return k_eff, rho_cp_eff
    
    def definir_condicoes_contorno(self):
        """
        Define condições de contorno dinâmicas
        """
        bcs = []
        
        # Condições de contorno para cada camada ativa
        for i in range(self.camada_atual + 1):
            camada_info = self.camadas_info[i]
            
            # Condições de Dirichlet (temperatura prescrita)
            if 'dirichlet' in camada_info:
                for bc_info in camada_info['dirichlet']:
                    facets = self.obter_facets_contorno(bc_info['tag'])
                    dofs = fem.locate_dofs_topological(self.V, self.mesh.topology.dim-1, facets)
                    bc_value = fem.Constant(self.mesh, bc_info['valor'])
                    bc = fem.dirichletbc(bc_value, dofs, self.V)
                    bcs.append(bc)
        
        return bcs
    
    def obter_facets_contorno(self, tag):
        """
        Retorna facetas do contorno com tag especificada
        """
        facets = []
        for i, facet_tag in enumerate(self.facet_tags.values):
            if facet_tag == tag:
                facets.append(i)
        return np.array(facets)
    
    def criar_formulacao_variacional(self, dt):
        """
        Cria formulação variacional para análise térmica transiente
        """
        # Função teste
        v = ufl.TestFunction(self.V)
        
        # Propriedades materiais efetivas
        k_eff, rho_cp_eff = self.definir_propriedades_material()
        
        # Formulação variacional
        # Termo transiente
        F = rho_cp_eff * (self.T - self.T_old) / dt * v * self.dx
        
        # Termo de difusão
        F += k_eff * ufl.dot(ufl.grad(self.T), ufl.grad(v)) * self.dx
        
        # Termo de geração interna
        F -= self.campo_geracao * v * self.dx
        
        # Condições de contorno de Neumann (fluxo prescrito)
        for i in range(self.camada_atual + 1):
            camada_info = self.camadas_info[i]
            if 'neumann' in camada_info:
                for bc_info in camada_info['neumann']:
                    tag = bc_info['tag']
                    fluxo = bc_info['fluxo']
                    F -= fluxo * v * self.ds(tag)
        
        # Condições de contorno de Robin (convecção)
        for i in range(self.camada_atual + 1):
            camada_info = self.camadas_info[i]
            if 'robin' in camada_info:
                for bc_info in camada_info['robin']:
                    tag = bc_info['tag']
                    h = bc_info['h']  # Coeficiente de convecção
                    T_amb = bc_info['T_amb']  # Temperatura ambiente
                    F += h * (self.T - T_amb) * v * self.ds(tag)
        
        return F
    
    def configurar_solver(self, dt):
        """
        Configura solver não-linear
        """
        # Formulação variacional
        F = self.criar_formulacao_variacional(dt)
        
        # Condições de contorno
        bcs = self.definir_condicoes_contorno()
        
        # Criar problema não-linear
        self.problem = fem.petsc.NonlinearProblem(F, self.T, bcs)
        
        # Configurar solver
        self.solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, self.problem)
        self.solver.convergence_criterion = "incremental"
        self.solver.rtol = 1e-6
        self.solver.atol = 1e-10
        self.solver.max_it = 25
        
        # Configurações do solver linear
        ksp = self.solver.krylov_solver
        opts = PETSc.Options()
        opts["ksp_type"] = "preonly"
        opts["pc_type"] = "lu"
        opts["pc_factor_mat_solver_type"] = "mumps"
        ksp.setFromOptions()
    
    def resolver_passo_tempo(self, dt):
        """
        Resolve um passo de tempo
        """
        # Atualizar solver se necessário
        if self.solver is None:
            self.configurar_solver(dt)
        
        # Resolver sistema não-linear
        n_iter, converged = self.solver.solve(self.T)
        
        if not converged:
            print(f"Aviso: Solver não convergiu em {n_iter} iterações")
            return False
        
        return True
    
    def simular_construcao(self, cronograma_construcao, dt_base=1.0):
        """
        Simula todo o processo de construção
        """
        print("Iniciando simulação de construção em camadas...")
        
        # Configurar arquivos de saída
        arquivo_saida = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "construcao_camadas.xdmf", "w")
        arquivo_saida.write_mesh(self.mesh)
        
        # Loop principal de simulação
        for evento in cronograma_construcao:
            tempo_alvo = evento['tempo']
            
            # Avançar no tempo até o próximo evento
            while self.tempo_atual < tempo_alvo:
                dt = min(dt_base, tempo_alvo - self.tempo_atual)
                
                # Resolver passo de tempo
                sucesso = self.resolver_passo_tempo(dt)
                if not sucesso:
                    print(f"Falha na convergência no tempo {self.tempo_atual}")
                    break
                
                # Atualizar solução anterior
                self.T_old.x.array[:] = self.T.x.array[:]
                self.tempo_atual += dt
                
                # Salvar resultado
                arquivo_saida.write_function(self.T, self.tempo_atual)
                
                print(f"Tempo: {self.tempo_atual:.2f}, T_max: {np.max(self.T.x.array):.2f}")
            
            # Processar evento
            if evento['tipo'] == 'nova_camada':
                nova_camada = evento['camada']
                self.atualizar_ativacao_camada(nova_camada)
                
                # Reconfigurar solver com nova formulação
                self.solver = None
                
                print(f"Camada {nova_camada} ativada no tempo {self.tempo_atual}")
        
        arquivo_saida.close()
        print("Simulação concluída!")
    
    def salvar_campos_ativacao(self, nome_arquivo):
        """
        Salva campos de ativação para debug
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, nome_arquivo, "w")
        arquivo.write_mesh(self.mesh)
        arquivo.write_function(self.campo_ativacao, 0.0)
        arquivo.write_function(self.campo_geracao, 0.0)
        arquivo.close()


def criar_exemplo_barragem():
    """
    Exemplo de uso para barragem
    """
    # Parâmetros da simulação
    parametros_termicos = {
        'condutividade': 2.5,      # W/m·K
        'densidade': 2400.0,       # kg/m³
        'calor_especifico': 1000.0 # J/kg·K
    }
    
    # Informações das camadas
    camadas_info = [
        {
            'tempo_lancamento': 0.0,
            'material': {'Q0': 1000.0, 'tau': 24.0},
            'dirichlet': [{'tag': 1, 'valor': 20.0}],
            'robin': [{'tag': 2, 'h': 10.0, 'T_amb': 15.0}]
        },
        {
            'tempo_lancamento': 168.0,  # 7 dias
            'material': {'Q0': 1000.0, 'tau': 24.0},
            'robin': [{'tag': 3, 'h': 10.0, 'T_amb': 15.0}]
        },
        {
            'tempo_lancamento': 336.0,  # 14 dias
            'material': {'Q0': 1000.0, 'tau': 24.0},
            'robin': [{'tag': 4, 'h': 10.0, 'T_amb': 15.0}]
        }
    ]
    
    # Cronograma de construção
    cronograma = [
        {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},
        {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},
        {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2},
        {'tempo': 720.0, 'tipo': 'fim_simulacao'}
    ]
    
    return parametros_termicos, camadas_info, cronograma


if __name__ == "__main__":
    # Exemplo de uso
    print("Exemplo de Construção em Camadas - FENICSx")
    
    # Criar malha (exemplo simples - na prática viria de arquivo .msh)
    mesh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, 50, 50)
    
    # Configurar simulação
    parametros, camadas, cronograma = criar_exemplo_barragem()
    
    # Criar objeto de simulação
    simulacao = ConstrucaoCamadas(mesh, camadas, parametros)
    
    # Executar simulação
    simulacao.simular_construcao(cronograma, dt_base=1.0)
    
    # Salvar campos de debug
    simulacao.salvar_campos_ativacao("campos_ativacao.xdmf")