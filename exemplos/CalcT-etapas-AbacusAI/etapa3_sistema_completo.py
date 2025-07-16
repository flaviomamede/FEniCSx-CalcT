
#!/usr/bin/python3
"""
ETAPA 3 - Sistema Completo
Implementa o sistema completo com:
- Espaços de funções para Tp, teq, Q
- Cálculo sequencial: teq -> Q -> Tp
- Funções update_equivalent_time() e calculate_heat_generation()
"""

import numpy as np
import json
from dolfinx import mesh, fem, io
from dolfinx.fem import LinearProblem
from mpi4py import MPI
import ufl
import petsc4py
petsc4py.init()
from petsc4py import PETSc

def update_equivalent_time(teq_anterior, dt, T_atual, T_ref=20.0, A=4000.0):
    """
    Atualiza o tempo equivalente usando a lei de Arrhenius
    teq_new = teq_old + dt * exp(A * (1/T_ref - 1/T_atual))
    """
    # Converter temperatura para Kelvin
    T_K = T_atual + 273.15
    T_ref_K = T_ref + 273.15
    
    # Calcular fator de Arrhenius
    arrhenius_factor = np.exp(A * (1.0/T_ref_K - 1.0/T_K))
    
    # Atualizar tempo equivalente
    teq_new = teq_anterior + dt * arrhenius_factor
    
    return teq_new

def calculate_heat_generation(teq, alpha_max=0.8, tau=86400.0, Q_max=50000.0):
    """
    Calcula a geração de calor baseada no tempo equivalente
    alpha = alpha_max * (1 - exp(-teq/tau))
    Q = Q_max * d(alpha)/dt = (Q_max * alpha_max / tau) * exp(-teq/tau)
    """
    # Calcular grau de hidratação
    alpha = alpha_max * (1.0 - np.exp(-teq / tau))
    
    # Calcular taxa de geração de calor
    Q = (Q_max * alpha_max / tau) * np.exp(-teq / tau)
    
    return Q, alpha

def main():
    print("=== ETAPA 3 - SISTEMA COMPLETO ===")
    
    # Parâmetros físicos
    rho = 2400.0  # kg/m³
    ce = 900.0    # J/(kg⋅K)
    k = 2.0       # W/(m⋅K)
    T_inicial = 20.0  # °C
    T_bc = 20.0       # °C
    dt = 3600.0       # s (1 hora)
    
    # Parâmetros de hidratação
    T_ref = 20.0      # °C (temperatura de referência)
    A = 4000.0        # K (parâmetro de Arrhenius)
    alpha_max = 0.8   # grau máximo de hidratação
    tau = 86400.0     # s (tempo característico - 24h)
    Q_max = 50000.0   # W/m³ (geração máxima de calor)
    
    print(f"Parâmetros físicos:")
    print(f"  ρ = {rho} kg/m³")
    print(f"  ce = {ce} J/(kg⋅K)")
    print(f"  k = {k} W/(m⋅K)")
    print(f"  T_inicial = {T_inicial} °C")
    print(f"  dt = {dt} s")
    
    print(f"Parâmetros de hidratação:")
    print(f"  T_ref = {T_ref} °C")
    print(f"  A = {A} K")
    print(f"  alpha_max = {alpha_max}")
    print(f"  tau = {tau} s")
    print(f"  Q_max = {Q_max} W/m³")
    
    # Carregar malha
    mesh_file = "/home/ubuntu/Uploads/barragem1.xdmf"
    print(f"\nCarregando malha: {mesh_file}")
    
    with io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
        domain = xdmf.read_mesh(name="malha")
        ct = xdmf.read_meshtags(domain, name="malha_cells")
    
    print(f"Malha carregada: {domain.topology.index_map(domain.topology.dim).size_global} células")
    
    # Carregar dados JSON
    json_file = "/home/ubuntu/Uploads/analise_stagewise_barragem1_xdmf.json"
    print(f"Carregando dados JSON: {json_file}")
    
    with open(json_file, 'r') as f:
        dados = json.load(f)
    
    # Extrair informações do Bloco 1
    bloco1 = dados["analise_resultados"]["bloco_1"]
    dominios_ativos = bloco1["physical_groups"]["surfaces"]
    nos_contorno = bloco1["elementos_nos"]["nos_contorno"]
    
    print(f"Bloco 1 - Domínios ativos: {dominios_ativos}")
    print(f"Bloco 1 - Nós de contorno: {nos_contorno}")
    
    # Definir espaços de funções
    V_T = fem.FunctionSpace(domain, ("Lagrange", 1))    # Temperatura
    V_teq = fem.FunctionSpace(domain, ("Lagrange", 1))  # Tempo equivalente
    V_Q = fem.FunctionSpace(domain, ("Lagrange", 1))    # Geração de calor
    
    print(f"Espaços de funções:")
    print(f"  V_T: {V_T.dofmap.index_map.size_global} DOFs")
    print(f"  V_teq: {V_teq.dofmap.index_map.size_global} DOFs")
    print(f"  V_Q: {V_Q.dofmap.index_map.size_global} DOFs")
    
    # Criar funções
    T_atual = fem.Function(V_T, name="T_atual")
    T_anterior = fem.Function(V_T, name="T_anterior")
    teq_atual = fem.Function(V_teq, name="teq_atual")
    teq_anterior = fem.Function(V_teq, name="teq_anterior")
    Q_atual = fem.Function(V_Q, name="Q_atual")
    
    # Inicializar campos
    T_anterior.x.array[:] = T_inicial
    T_atual.x.array[:] = T_inicial
    teq_anterior.x.array[:] = 0.0  # Tempo equivalente inicial = 0
    teq_atual.x.array[:] = 0.0
    Q_atual.x.array[:] = 0.0
    
    print(f"Campos inicializados:")
    print(f"  T_inicial = {T_inicial} °C")
    print(f"  teq_inicial = 0.0 s")
    print(f"  Q_inicial = 0.0 W/m³")
    
    # Definir condições de contorno para temperatura
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
    )
    
    boundary_dofs = fem.locate_dofs_topological(V_T, fdim, boundary_facets)
    
    print(f"DOFs de contorno encontrados: {len(boundary_dofs)}")
    
    # Criar BC de Dirichlet para temperatura
    T_bc_func = fem.Function(V_T)
    T_bc_func.x.array[:] = T_bc
    
    bc_T = fem.DirichletBC(T_bc_func, boundary_dofs)
    bcs = [bc_T]
    
    print(f"Condição de contorno Dirichlet aplicada: T = {T_bc} °C")
    
    # Definir funções de teste e trial para temperatura
    u = ufl.TrialFunction(V_T)
    v = ufl.TestFunction(V_T)
    
    # Medidas de integração
    dx = ufl.Measure("dx", domain=domain, subdomain_data=ct)
    
    print("\nDefinindo forma variacional para temperatura...")
    
    # Forma bilinear para temperatura
    a = (rho * ce / dt) * u * v * ufl.dx
    a += k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    
    print("Forma variacional definida")
    
    # Resolver sistema acoplado
    print("\n=== INICIANDO SIMULAÇÃO ACOPLADA ===")
    
    # Executar 5 passos de tempo
    for step in range(1, 6):
        print(f"\nPasso {step}:")
        print(f"  Tempo = {step * dt} s")
        
        # ETAPA 3.1: Atualizar tempo equivalente
        print("  3.1 - Atualizando tempo equivalente...")
        teq_new = update_equivalent_time(teq_anterior.x.array, dt, T_anterior.x.array, T_ref, A)
        teq_atual.x.array[:] = teq_new
        
        teq_min = np.min(teq_atual.x.array)
        teq_max = np.max(teq_atual.x.array)
        teq_mean = np.mean(teq_atual.x.array)
        print(f"    teq_min = {teq_min:.2f} s")
        print(f"    teq_max = {teq_max:.2f} s")
        print(f"    teq_mean = {teq_mean:.2f} s")
        
        # ETAPA 3.2: Calcular geração de calor
        print("  3.2 - Calculando geração de calor...")
        Q_new, alpha = calculate_heat_generation(teq_atual.x.array, alpha_max, tau, Q_max)
        Q_atual.x.array[:] = Q_new
        
        Q_min = np.min(Q_atual.x.array)
        Q_max_val = np.max(Q_atual.x.array)
        Q_mean = np.mean(Q_atual.x.array)
        alpha_mean = np.mean(alpha)
        print(f"    Q_min = {Q_min:.2f} W/m³")
        print(f"    Q_max = {Q_max_val:.2f} W/m³")
        print(f"    Q_mean = {Q_mean:.2f} W/m³")
        print(f"    alpha_mean = {alpha_mean:.4f}")
        
        # ETAPA 3.3: Resolver equação de temperatura
        print("  3.3 - Resolvendo equação de temperatura...")
        
        # Forma linear com geração de calor atual
        L = (rho * ce / dt) * T_anterior * v * ufl.dx
        L += Q_atual * v * ufl.dx  # Adicionar geração de calor
        
        # Criar e resolver problema linear
        problem = LinearProblem(a, L, bcs=bcs, 
                              petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        
        try:
            T_atual = problem.solve()
            print(f"    Solução obtida com sucesso")
            
            # Diagnósticos de temperatura
            T_min = np.min(T_atual.x.array)
            T_max = np.max(T_atual.x.array)
            T_mean = np.mean(T_atual.x.array)
            
            print(f"    T_min = {T_min:.6f} °C")
            print(f"    T_max = {T_max:.6f} °C")
            print(f"    T_mean = {T_mean:.6f} °C")
            
            # Verificar se solução é finita
            if np.all(np.isfinite(T_atual.x.array)):
                print(f"    ✓ Solução é finita")
            else:
                print(f"    ✗ Solução contém valores não-finitos!")
                return False
            
            # Verificar se não divergiu para infinito
            if T_max < 1000.0:  # Limite razoável
                print(f"    ✓ Temperatura dentro de limites razoáveis")
            else:
                print(f"    ✗ Temperatura muito alta - possível divergência!")
                return False
            
            # Atualizar campos anteriores
            T_anterior.x.array[:] = T_atual.x.array[:]
            teq_anterior.x.array[:] = teq_atual.x.array[:]
            
        except Exception as e:
            print(f"    ✗ Erro na solução: {e}")
            return False
    
    # Salvar resultados
    print(f"\nSalvando resultados...")
    
    try:
        # Salvar temperatura
        with io.XDMFFile(MPI.COMM_WORLD, "/home/ubuntu/etapa3_temperatura.xdmf", "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(T_atual)
        
        # Salvar tempo equivalente
        with io.XDMFFile(MPI.COMM_WORLD, "/home/ubuntu/etapa3_teq.xdmf", "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(teq_atual)
        
        # Salvar geração de calor
        with io.XDMFFile(MPI.COMM_WORLD, "/home/ubuntu/etapa3_Q.xdmf", "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(Q_atual)
        
        print("✓ Resultados salvos com sucesso")
    except Exception as e:
        print(f"✗ Erro ao salvar: {e}")
    
    print("\n=== ETAPA 3 CONCLUÍDA COM SUCESSO ===")
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("RESULTADO: SUCESSO")
        exit(0)
    else:
        print("RESULTADO: FALHA")
        exit(1)
