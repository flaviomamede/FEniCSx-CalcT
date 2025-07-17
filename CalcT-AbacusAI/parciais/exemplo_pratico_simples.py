#!/usr/bin/env python3
"""
Exemplo Prático Simples - Construção em Camadas FENICSx
Versão simplificada para facilitar o início da implementação
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.io
import dolfinx.nls.petsc
import ufl
from mpi4py import MPI
import time

class ConstrucaoCamadasSimples:
    """
    Versão simplificada da construção em camadas
    """
    
    def __init__(self, mesh, parametros_termicos):
        self.mesh = mesh
        self.parametros = parametros_termicos
        
        # Espaço funcional
        self.V = fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solução
        self.T = fem.Function(self.V)
        self.T_old = fem.Function(self.V)
        
        # Campos de controle
        self.ativacao = fem.Function(self.V)
        self.geracao = fem.Function(self.V)
        
        # Medidas
        self.dx = ufl.Measure("dx", domain=mesh)
        self.ds = ufl.Measure("ds", domain=mesh)
        
        # Estado
        self.tempo = 0.0
        self.camadas_ativas = set()
        
        # Inicializar
        self.T.x.array[:] = 20.0  # Temperatura inicial
        self.T_old.x.array[:] = 20.0
        self.ativacao.x.array[:] = 0.0
        self.geracao.x.array[:] = 0.0
        
        print("Simulador simples inicializado com sucesso!")
        print(f"DOFs: {self.V.dofmap.index_map.size_local}")
    
    def ativar_camada(self, num_camada, Q0, tau):
        """
        Ativa uma camada específica
        """
        if num_camada in self.camadas_ativas:
            return
        
        print(f"Ativando camada {num_camada}")
        
        # Definir região da camada (exemplo: divisão vertical)
        height = 1.0
        layer_height = height / 3  # 3 camadas
        y_min = num_camada * layer_height
        y_max = (num_camada + 1) * layer_height
        
        # Ativar DOFs na região da camada
        for i in range(self.V.dofmap.index_map.size_local):
            coord = self.V.tabulate_dof_coordinates()[i]
            if y_min <= coord[1] <= y_max:
                self.ativacao.x.array[i] = 1.0
                self.geracao.x.array[i] = Q0 * np.exp(-self.tempo / tau)
        
        self.camadas_ativas.add(num_camada)
        
        # Atualizar campos
        self.ativacao.x.scatter_forward()
        self.geracao.x.scatter_forward()
        
        print(f"Camada {num_camada} ativada com {np.sum(self.ativacao.x.array > 0.5)} DOFs")
    
    def atualizar_geracao(self, parametros_camadas):
        """
        Atualiza geração de calor com base no tempo
        """
        for num_camada in self.camadas_ativas:
            if num_camada in parametros_camadas:
                Q0 = parametros_camadas[num_camada]['Q0']
                tau = parametros_camadas[num_camada]['tau']
                tempo_lancamento = parametros_camadas[num_camada]['tempo_lancamento']
                
                # Calcular idade da camada
                idade = max(0.0, self.tempo - tempo_lancamento)
                if idade > 0:
                    geracao_atual = Q0 * np.exp(-idade / tau)
                    
                    # Atualizar apenas DOFs desta camada
                    height = 1.0
                    layer_height = height / 3
                    y_min = num_camada * layer_height
                    y_max = (num_camada + 1) * layer_height
                    
                    for i in range(self.V.dofmap.index_map.size_local):
                        coord = self.V.tabulate_dof_coordinates()[i]
                        if y_min <= coord[1] <= y_max and self.ativacao.x.array[i] > 0.5:
                            self.geracao.x.array[i] = geracao_atual
        
        self.geracao.x.scatter_forward()
    
    def criar_formulacao(self, dt):
        """
        Cria formulação variacional
        """
        v = ufl.TestFunction(self.V)
        
        # Propriedades efetivas
        k = self.ativacao * self.parametros['k']
        rho_cp = self.ativacao * self.parametros['rho'] * self.parametros['cp']
        
        # Formulação variacional
        F = (rho_cp * (self.T - self.T_old) / dt) * v * self.dx
        F += k * ufl.dot(ufl.grad(self.T), ufl.grad(v)) * self.dx
        F -= self.geracao * v * self.dx
        
        return F
    
    def definir_contorno(self):
        """
        Define condições de contorno simples
        """
        # Contorno inferior: temperatura fixa
        def contorno_inferior(x):
            return np.isclose(x[1], 0.0)
        
        # Localizar DOFs no contorno
        dofs_contorno = fem.locate_dofs_geometrical(self.V, contorno_inferior)
        
        # Condição de Dirichlet
        bc = fem.dirichletbc(fem.Constant(self.mesh, 20.0), dofs_contorno, self.V)
        
        return [bc]
    
    def resolver_passo(self, dt):
        """
        Resolve um passo de tempo
        """
        # Criar formulação
        F = self.criar_formulacao(dt)
        
        # Condições de contorno
        bcs = self.definir_contorno()
        
        # Criar problema
        problema = fem.petsc.NonlinearProblem(F, self.T, bcs)
        
        # Solver
        solver = dolfinx.nls.petsc.NewtonSolver(MPI.COMM_WORLD, problema)
        solver.rtol = 1e-6
        solver.atol = 1e-10
        solver.max_it = 25
        
        # Resolver
        n_iter, converged = solver.solve(self.T)
        
        # Atualizar solução anterior
        self.T_old.x.array[:] = self.T.x.array.copy()
        
        return converged
    
    def simular(self, cronograma, dt=2.0):
        """
        Executa simulação completa
        """
        print("="*50)
        print("SIMULAÇÃO INICIADA")
        print("="*50)
        
        # Arquivo de saída
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "simulacao_simples.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        
        # Configurar nomes
        self.T.name = "Temperatura"
        self.ativacao.name = "Ativacao"
        self.geracao.name = "Geracao"
        
        # Parâmetros das camadas
        parametros_camadas = {}
        
        # Processar cronograma
        for evento in cronograma:
            tempo_alvo = evento['tempo']
            
            print(f"\nProcessando até t={tempo_alvo:.1f}h")
            
            # Avançar no tempo
            while self.tempo < tempo_alvo:
                dt_atual = min(dt, tempo_alvo - self.tempo)
                
                # Atualizar geração
                self.atualizar_geracao(parametros_camadas)
                
                # Resolver passo
                converged = self.resolver_passo(dt_atual)
                
                if not converged:
                    print(f"Aviso: Não convergiu no tempo {self.tempo:.1f}h")
                
                # Atualizar tempo
                self.tempo += dt_atual
                
                # Salvar resultado
                arquivo.write_function(self.T, self.tempo)
            
            # Processar evento
            if evento['tipo'] == 'camada':
                parametros_camadas[evento['num']] = {
                    'Q0': evento['Q0'],
                    'tau': evento['tau'],
                    'tempo_lancamento': self.tempo
                }
                self.ativar_camada(evento['num'], evento['Q0'], evento['tau'])
            
            elif evento['tipo'] == 'fim':
                print("Simulação concluída!")
                break
        
        arquivo.close()
        
        # Salvar campos finais
        self.salvar_campos_finais()
        
        print(f"\nSimulação concluída!")
        print(f"Tempo total: {self.tempo:.1f}h")
        print(f"Camadas ativas: {sorted(self.camadas_ativas)}")
        
        # Estatísticas finais
        T_min = np.min(self.T.x.array)
        T_max = np.max(self.T.x.array)
        T_med = np.mean(self.T.x.array)
        
        print(f"Temperatura final: {T_min:.1f}°C - {T_max:.1f}°C (média: {T_med:.1f}°C)")
        
        return self
    
    def salvar_campos_finais(self):
        """
        Salva campos finais para análise
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "campos_finais_simples.xdmf", "w")
        arquivo.write_mesh(self.mesh)
        arquivo.write_function(self.T, self.tempo)
        arquivo.write_function(self.ativacao, self.tempo)
        arquivo.write_function(self.geracao, self.tempo)
        arquivo.close()
        
        print("Campos finais salvos em: campos_finais_simples.xdmf")


def exemplo_uso():
    """
    Exemplo de uso básico
    """
    print("Exemplo Simplificado - Construção em Camadas")
    print("="*50)
    
    # Criar malha retangular simples
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [2.0, 1.0]],
        [20, 10]
    )
    
    print(f"Malha criada: {mesh.topology.index_map(2).size_local} células")
    
    # Parâmetros térmicos
    parametros = {
        'k': 2.5,      # Condutividade térmica [W/m·K]
        'rho': 2400.0, # Densidade [kg/m³]
        'cp': 1000.0   # Calor específico [J/kg·K]
    }
    
    # Criar simulador
    sim = ConstrucaoCamadasSimples(mesh, parametros)
    
    # Cronograma simples
    cronograma = [
        {'tempo': 0.0, 'tipo': 'camada', 'num': 0, 'Q0': 1000.0, 'tau': 48.0},
        {'tempo': 24.0, 'tipo': 'camada', 'num': 1, 'Q0': 800.0, 'tau': 48.0},
        {'tempo': 48.0, 'tipo': 'camada', 'num': 2, 'Q0': 600.0, 'tau': 48.0},
        {'tempo': 168.0, 'tipo': 'fim'}
    ]
    
    # Executar simulação
    sim.simular(cronograma, dt=2.0)
    
    print("\nArquivos gerados:")
    print("- simulacao_simples.xdmf (evolução temporal)")
    print("- campos_finais_simples.xdmf (campos finais)")
    print("\nUse ParaView para visualizar os resultados!")
    
    return sim


if __name__ == "__main__":
    exemplo_uso()