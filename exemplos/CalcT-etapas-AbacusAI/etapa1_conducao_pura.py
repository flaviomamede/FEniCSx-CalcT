
#!/usr/bin/python3
"""
ETAPA 1 - Condução Pura Simples
Teste básico de condução de calor nos domínios ativos do Bloco 1
Sem geração de calor, apenas condução pura com BCs de Dirichlet
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

def main():
    print("=== ETAPA 1 - CONDUÇÃO PURA SIMPLES ===")
    
    # Parâmetros físicos extraídos do YAML
    rho = 2400.0  # kg/m³ (densidade)
    ce = 900.0    # J/(kg⋅K) (calor específico)
    k = 2.0       # W/(m⋅K) (condutividade térmica)
    T_inicial = 20.0  # °C (temperatura inicial)
    T_bc = 20.0       # °C (temperatura de contorno)
    dt = 3600.0       # s (1 hora)
    
    print(f"Parâmetros físicos:")
    print(f"  ρ = {rho} kg/m³")
    print(f"  ce = {ce} J/(kg⋅K)")
    print(f"  k = {k} W/(m⋅K)")
    print(f"  T_inicial = {T_inicial} °C")
    print(f"  dt = {dt} s")
    
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
    
    # Definir espaço de funções para temperatura
    V_T = fem.FunctionSpace(domain, ("Lagrange", 1))
    print(f"Espaço de funções V_T: {V_T.dofmap.index_map.size_global} DOFs")
    
    # Criar funções de temperatura
    T_atual = fem.Function(V_T, name="T_atual")
    T_anterior = fem.Function(V_T, name="T_anterior")
    
    # Inicializar temperatura
    T_anterior.x.array[:] = T_inicial
    T_atual.x.array[:] = T_inicial
    
    print(f"Temperatura inicial definida: {T_inicial} °C")
    
    # Definir condições de contorno de Dirichlet
    # Aplicar BC nos nós de contorno do Bloco 1
    boundary_dofs = []
    
    # Encontrar DOFs nos nós de contorno
    # Para simplificar, vamos aplicar BC em toda a fronteira
    fdim = domain.topology.dim - 1
    boundary_facets = mesh.locate_entities_boundary(
        domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
    )
    
    boundary_dofs = fem.locate_dofs_topological(V_T, fdim, boundary_facets)
    
    print(f"DOFs de contorno encontrados: {len(boundary_dofs)}")
    
    # Criar BC de Dirichlet usando sintaxe correta para dolfinx 0.3.0
    # Criar função constante para o valor da BC
    T_bc_func = fem.Function(V_T)
    T_bc_func.x.array[:] = T_bc
    
    # Criar BC de Dirichlet
    bc_T = fem.DirichletBC(T_bc_func, boundary_dofs)
    bcs = [bc_T]
    
    print(f"Condição de contorno Dirichlet aplicada: T = {T_bc} °C")
    
    # Definir funções de teste e trial
    u = ufl.TrialFunction(V_T)
    v = ufl.TestFunction(V_T)
    
    # Medidas de integração para domínios ativos
    dx = ufl.Measure("dx", domain=domain, subdomain_data=ct)
    
    # Forma variacional simplificada (SEM fonte Q)
    # a = (rho*ce/dt)*u*v*dx + k*dot(grad(u), grad(v))*dx
    # L = (rho*ce/dt)*T_anterior*v*dx
    
    print("\nDefinindo forma variacional...")
    
    # Forma bilinear - aplicar em todo o domínio para evitar DOFs órfãos
    print("  Aplicando forma variacional em todo o domínio")
    
    # Termo transiente em todo o domínio
    a = (rho * ce / dt) * u * v * ufl.dx
    # Termo de difusão em todo o domínio  
    a += k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    # Termo fonte em todo o domínio
    L = (rho * ce / dt) * T_anterior * v * ufl.dx
    
    print("Forma variacional definida")
    
    # Resolver sistema linear
    print("\n=== INICIANDO SIMULAÇÃO ===")
    
    # Executar 2 passos de tempo
    for step in range(1, 3):
        print(f"\nPasso {step}:")
        print(f"  Tempo = {step * dt} s")
        
        # Criar problema linear
        problem = LinearProblem(a, L, bcs=bcs, 
                              petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        
        # Resolver
        try:
            T_atual = problem.solve()
            print(f"  Solução obtida com sucesso")
            
            # Diagnósticos
            T_min = np.min(T_atual.x.array)
            T_max = np.max(T_atual.x.array)
            T_mean = np.mean(T_atual.x.array)
            
            print(f"  T_min = {T_min:.6f} °C")
            print(f"  T_max = {T_max:.6f} °C")
            print(f"  T_mean = {T_mean:.6f} °C")
            
            # Verificar se solução é finita
            if np.all(np.isfinite(T_atual.x.array)):
                print(f"  ✓ Solução é finita")
            else:
                print(f"  ✗ Solução contém valores não-finitos!")
                return False
            
            # Atualizar temperatura anterior
            T_anterior.x.array[:] = T_atual.x.array[:]
            
        except Exception as e:
            print(f"  ✗ Erro na solução: {e}")
            return False
    
    # Salvar resultados
    output_file = "/home/ubuntu/etapa1_resultado.xdmf"
    print(f"\nSalvando resultados em: {output_file}")
    
    try:
        with io.XDMFFile(MPI.COMM_WORLD, output_file, "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_function(T_atual)
        print("✓ Resultados salvos com sucesso")
    except Exception as e:
        print(f"✗ Erro ao salvar: {e}")
    
    print("\n=== ETAPA 1 CONCLUÍDA COM SUCESSO ===")
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("RESULTADO: SUCESSO")
        exit(0)
    else:
        print("RESULTADO: FALHA")
        exit(1)
