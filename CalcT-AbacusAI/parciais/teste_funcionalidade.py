#!/usr/bin/env python3
"""
Exemplo Teste Básico - Construção em Camadas FENICSx
Teste simples para verificar funcionamento
"""

import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.fem.petsc
import dolfinx.io
import ufl
from mpi4py import MPI


def teste_basico():
    """
    Teste básico da construção em camadas
    """
    print("="*60)
    print("TESTE BÁSICO - CONSTRUÇÃO EM CAMADAS")
    print("="*60)
    
    # Criar malha simples
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [4.0, 3.0]],
        [16, 12]
    )
    
    print(f"Malha: {mesh.topology.index_map(2).size_local} células")
    
    # Espaço funcional
    V = fem.functionspace(mesh, ("Lagrange", 1))
    print(f"DOFs: {V.dofmap.index_map.size_local}")
    
    # Campos
    T = fem.Function(V)
    T_old = fem.Function(V)
    ativacao = fem.Function(V)
    geracao = fem.Function(V)
    
    # Inicializar
    T.x.array[:] = 20.0
    T_old.x.array[:] = 20.0
    ativacao.x.array[:] = 0.0
    geracao.x.array[:] = 0.0
    
    print("Campos inicializados")
    
    # Parâmetros
    k = 2.5
    rho = 2400.0
    cp = 1000.0
    dt = 3600.0  # 1 hora
    
    # Medidas
    dx = ufl.Measure("dx", domain=mesh)
    
    # Ativar primeira camada (terço inferior)
    coords = V.tabulate_dof_coordinates()
    for i in range(len(coords)):
        if coords[i][1] < 1.0:  # Primeira camada
            ativacao.x.array[i] = 1.0
            geracao.x.array[i] = 500.0  # Geração de calor
    
    ativacao.x.scatter_forward()
    geracao.x.scatter_forward()
    
    print(f"Primeira camada ativada: {np.sum(ativacao.x.array > 0.5)} DOFs")
    
    # Condições de contorno
    def base(x):
        return np.isclose(x[1], 0.0)
    
    dofs_base = fem.locate_dofs_geometrical(V, base)
    bc = fem.dirichletbc(fem.Constant(mesh, 20.0), dofs_base, V)
    
    print(f"Condição de contorno: {len(dofs_base)} DOFs")
    
    # Formulação linear simples
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    # Forma bilinear (estacionária)
    a = (ativacao * k * ufl.dot(ufl.grad(u), ufl.grad(v))) * dx
    
    # Forma linear
    L = geracao * v * dx
    
    print("Formulação criada")
    
    # Resolver problema linear
    problema = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc])
    T_resultado = problema.solve()
    
    print("Problema resolvido!")
    
    # Estatísticas
    T_min = np.min(T_resultado.x.array)
    T_max = np.max(T_resultado.x.array)
    T_med = np.mean(T_resultado.x.array)
    
    print(f"Temperatura: {T_min:.1f}°C - {T_max:.1f}°C (média: {T_med:.1f}°C)")
    
    # Salvar resultado
    T_resultado.name = "Temperatura"
    ativacao.name = "Ativacao"
    geracao.name = "Geracao"
    
    arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "teste_basico.xdmf", "w")
    arquivo.write_mesh(mesh)
    arquivo.write_function(T_resultado, 0.0)
    arquivo.write_function(ativacao, 0.0)
    arquivo.write_function(geracao, 0.0)
    arquivo.close()
    
    print("Resultado salvo em: teste_basico.xdmf")
    
    return T_resultado, ativacao, geracao


def teste_transiente():
    """
    Teste transiente simples
    """
    print("\n" + "="*60)
    print("TESTE TRANSIENTE")
    print("="*60)
    
    # Criar malha
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [2.0, 1.0]],
        [10, 5]
    )
    
    # Espaço funcional
    V = fem.functionspace(mesh, ("Lagrange", 1))
    
    # Campos
    T = fem.Function(V)
    T_old = fem.Function(V)
    
    # Inicializar
    T.x.array[:] = 20.0
    T_old.x.array[:] = 20.0
    
    # Parâmetros
    k = 2.5
    rho = 2400.0
    cp = 1000.0
    dt = 3600.0  # 1 hora
    Q = 1000.0   # Geração constante
    
    # Condições de contorno
    def contorno(x):
        return np.isclose(x[1], 0.0)
    
    dofs_contorno = fem.locate_dofs_geometrical(V, contorno)
    bc = fem.dirichletbc(fem.Constant(mesh, 20.0), dofs_contorno, V)
    
    # Formulação transiente
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    dx = ufl.Measure("dx", domain=mesh)
    
    # Forma bilinear
    a = (rho * cp / dt) * u * v * dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx
    
    # Forma linear
    L = (rho * cp / dt) * T_old * v * dx + Q * v * dx
    
    # Arquivo de saída
    arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "teste_transiente.xdmf", "w")
    arquivo.write_mesh(mesh)
    T.name = "Temperatura"
    
    # Loop de tempo
    tempo = 0.0
    tempo_final = 24.0  # 24 horas
    
    print(f"Simulação transiente: 0 - {tempo_final} horas")
    
    while tempo < tempo_final:
        # Resolver
        problema = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc])
        T_new = problema.solve()
        
        # Atualizar
        T.x.array[:] = T_new.x.array.copy()
        T_old.x.array[:] = T_new.x.array.copy()
        
        tempo += dt / 3600.0  # Converter para horas
        
        # Salvar
        arquivo.write_function(T, tempo)
        
        # Status
        if int(tempo) % 6 == 0:  # A cada 6 horas
            T_med = np.mean(T.x.array)
            print(f"t={tempo:4.1f}h: T_média={T_med:.1f}°C")
    
    arquivo.close()
    
    print("Teste transiente concluído!")
    print("Resultado salvo em: teste_transiente.xdmf")
    
    return T


def teste_camadas():
    """
    Teste com ativação de camadas
    """
    print("\n" + "="*60)
    print("TESTE COM CAMADAS")
    print("="*60)
    
    # Criar malha
    mesh = dolfinx.mesh.create_rectangle(
        MPI.COMM_WORLD,
        [[0.0, 0.0], [4.0, 3.0]],
        [20, 15]
    )
    
    # Espaço funcional
    V = fem.functionspace(mesh, ("Lagrange", 1))
    
    # Campos
    T = fem.Function(V)
    T_old = fem.Function(V)
    ativacao = fem.Function(V)
    geracao = fem.Function(V)
    
    # Inicializar
    T.x.array[:] = 20.0
    T_old.x.array[:] = 20.0
    ativacao.x.array[:] = 0.0
    geracao.x.array[:] = 0.0
    
    # Parâmetros
    k = 2.5
    rho = 2400.0
    cp = 1000.0
    dt = 3600.0  # 1 hora
    
    # Condições de contorno
    def base(x):
        return np.isclose(x[1], 0.0)
    
    dofs_base = fem.locate_dofs_geometrical(V, base)
    bc = fem.dirichletbc(fem.Constant(mesh, 20.0), dofs_base, V)
    
    # Arquivo de saída
    arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, "teste_camadas.xdmf", "w")
    arquivo.write_mesh(mesh)
    T.name = "Temperatura"
    ativacao.name = "Ativacao"
    geracao.name = "Geracao"
    
    # Coordenadas dos DOFs
    coords = V.tabulate_dof_coordinates()
    
    # Simulação com 3 camadas
    tempo = 0.0
    camadas_ativadas = set()
    
    cronograma = [
        {'tempo': 0.0, 'camada': 0, 'y_min': 0.0, 'y_max': 1.0, 'Q': 1000.0},
        {'tempo': 24.0, 'camada': 1, 'y_min': 1.0, 'y_max': 2.0, 'Q': 800.0},
        {'tempo': 48.0, 'camada': 2, 'y_min': 2.0, 'y_max': 3.0, 'Q': 600.0}
    ]
    
    for evento in cronograma:
        tempo_alvo = evento['tempo']
        
        print(f"\nAtivar camada {evento['camada']} no tempo {tempo_alvo}h")
        
        # Ativar camada
        if evento['camada'] not in camadas_ativadas:
            for i in range(len(coords)):
                if evento['y_min'] <= coords[i][1] < evento['y_max']:
                    ativacao.x.array[i] = 1.0
                    geracao.x.array[i] = evento['Q']
            
            camadas_ativadas.add(evento['camada'])
            ativacao.x.scatter_forward()
            geracao.x.scatter_forward()
        
        # Simular até tempo alvo
        while tempo < tempo_alvo + 24.0:  # +24h após ativação
            # Formulação
            u = ufl.TrialFunction(V)
            v = ufl.TestFunction(V)
            dx = ufl.Measure("dx", domain=mesh)
            
            # Propriedades efetivas
            k_eff = ativacao * k + (1 - ativacao) * 0.01  # Pequena condutividade para elementos inativos
            rho_cp_eff = ativacao * rho * cp + (1 - ativacao) * 1.0  # Pequena capacidade
            
            # Forma bilinear
            a = (rho_cp_eff / dt) * u * v * dx + k_eff * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx
            
            # Forma linear
            L = (rho_cp_eff / dt) * T_old * v * dx + geracao * v * dx
            
            # Resolver
            problema = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc])
            T_new = problema.solve()
            
            # Atualizar
            T.x.array[:] = T_new.x.array.copy()
            T_old.x.array[:] = T_new.x.array.copy()
            
            tempo += dt / 3600.0
            
            # Salvar
            arquivo.write_function(T, tempo)
            arquivo.write_function(ativacao, tempo)
            arquivo.write_function(geracao, tempo)
            
            # Status
            if int(tempo) % 12 == 0:
                T_med = np.mean(T.x.array)
                dofs_ativos = np.sum(ativacao.x.array > 0.5)
                print(f"t={tempo:4.1f}h: T_média={T_med:.1f}°C, DOFs ativos={dofs_ativos}")
    
    arquivo.close()
    
    print(f"\nTeste com camadas concluído!")
    print(f"Camadas ativadas: {sorted(camadas_ativadas)}")
    print("Resultado salvo em: teste_camadas.xdmf")
    
    return T


if __name__ == "__main__":
    # Executar testes
    print("Executando testes de funcionalidade...")
    
    # Teste 1: Básico
    T1, ativacao1, geracao1 = teste_basico()
    
    # Teste 2: Transiente
    T2 = teste_transiente()
    
    # Teste 3: Camadas
    T3 = teste_camadas()
    
    print("\n" + "="*60)
    print("TODOS OS TESTES EXECUTADOS COM SUCESSO!")
    print("="*60)
    print("Arquivos gerados:")
    print("- teste_basico.xdmf")
    print("- teste_transiente.xdmf")
    print("- teste_camadas.xdmf")
    print("\nUse ParaView para visualizar os resultados!")