#!/usr/bin/python3
"""
ETAPA 4 - Sistema com Submeshes
Baseado na ETAPA 3, mas preparado para trabalhar com submeshes:
- Espaços de funções para Tp, teq, Q
- Cálculo sequencial: teq -> Q -> Tp
- Preparado para implementação de submeshes por etapas construtivas
- Funções update_equivalent_time() e calculate_heat_generation()
"""

import numpy as np
import json
import os
from dolfinx import mesh, fem, io
from mpi4py import MPI
import ufl
import petsc4py
petsc4py.init()
from petsc4py import PETSc

def update_equivalent_time(teq_anterior, dt, T_atual, T_ref=20.0, A=4000.0):
    """
    Atualiza o tempo equivalente usando a lei de Arrhenius
    teq_new = teq_old + dt * exp(A * (1/T_ref - 1/T_atual))
    
    Parâmetros:
    - teq_anterior: tempo equivalente do passo anterior (s)
    - dt: passo de tempo (s)
    - T_atual: temperatura atual (°C)
    - T_ref: temperatura de referência (°C)
    - A: parâmetro de Arrhenius (K)
    """
    # Converter temperatura para Kelvin
    T_K = T_atual + 273.15
    T_ref_K = T_ref + 273.15
    
    # Calcular fator de Arrhenius
    arrhenius_factor = np.exp(A * (1.0/T_ref_K - 1.0/T_K))
    
    # Atualizar tempo equivalente
    teq_new = teq_anterior + dt * arrhenius_factor
    
    return teq_new

def calculate_heat_generation(teq, alpha_max=0.8, tau=86400.0, Q_total=300000000.0):
    """
    Calcula a geração de calor baseada no tempo equivalente
    
    Modelo de hidratação:
    alpha = alpha_max * (1 - exp(-teq/tau))
    Q = Q_total * d(alpha)/dt = (Q_total * alpha_max / tau) * exp(-teq/tau)
    
    Parâmetros:
    - teq: tempo equivalente (s)
    - alpha_max: grau máximo de hidratação (adimensional)
    - tau: tempo característico (s)
    - Q_total: calor total de hidratação (J/m³)
    """
    # Calcular grau de hidratação
    alpha = alpha_max * (1.0 - np.exp(-teq / tau))
    
    # Calcular taxa de geração de calor
    Q = (Q_total * alpha_max / tau) * np.exp(-teq / tau)
    
    return Q, alpha

def main():
    print("=== ETAPA 4 - SISTEMA COM SUBMESHES ===")
    
    # Parâmetros físicos
    rho = 2400.0  # kg/m³ - densidade do concreto
    ce = 900.0    # J/(kg⋅K) - calor específico
    k = 2.0       # W/(m⋅K) - condutividade térmica
    T_inicial = 20.0  # °C - temperatura inicial
    T_bc = 20.0       # °C - temperatura de contorno
    
    # Parâmetros de hidratação
    T_ref = 20.0      # °C (temperatura de referência)
    A = 4000.0        # K (parâmetro de Arrhenius)
    alpha_max = 0.8   # grau máximo de hidratação
    tau = 86400.0     # s (tempo característico - 24h)
    Q_total = 300000000.0  # J/m³ (calor total de hidratação - 300 MJ/m³)
    
    print(f"Parâmetros físicos:")
    print(f"  ρ = {rho} kg/m³")
    print(f"  ce = {ce} J/(kg⋅K)")
    print(f"  k = {k} W/(m⋅K)")
    print(f"  T_inicial = {T_inicial} °C")
    
    print(f"Parâmetros de hidratação:")
    print(f"  T_ref = {T_ref} °C")
    print(f"  A = {A} K")
    print(f"  alpha_max = {alpha_max}")
    print(f"  tau = {tau} s")
    print(f"  Q_total = {Q_total/1e6:.1f} MJ/m³")
    
    # Caminhos dos arquivos - usando malha refinada
    base_path = "Uploads"
    mesh_file = os.path.join(base_path, "barragem2.xdmf")
    json_file = os.path.join(base_path, "analise_stagewise_barragem2_xdmf.json")
    
    # Verificar se arquivos existem
    if not os.path.exists(mesh_file):
        print(f"ERRO: Arquivo de malha não encontrado: {mesh_file}")
        return False
    
    if not os.path.exists(json_file):
        print(f"ERRO: Arquivo JSON não encontrado: {json_file}")
        return False
    
    # Carregar malha principal
    print(f"\nCarregando malha principal: {mesh_file}")
    
    try:
        with io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
            mesh_global = xdmf.read_mesh(name="malha")
            ct_global = xdmf.read_meshtags(mesh_global, name="malha_cells")
        
        print(f"Malha global carregada: {mesh_global.topology.index_map(mesh_global.topology.dim).size_global} células")
        
    except Exception as e:
        print(f"ERRO ao carregar malha: {e}")
        return False
    
    # Carregar dados JSON
    print(f"Carregando dados JSON: {json_file}")
    
    try:
        with open(json_file, 'r') as f:
            dados = json.load(f)
        
        # Extrair informações do Bloco 1
        bloco1 = dados["analise_resultados"]["bloco_1"]
        dominios_ativos = bloco1["physical_groups"]["surfaces"]
        nos_contorno = bloco1["elementos_nos"]["nos_contorno"]
        elementos_ativos_bloco1 = bloco1["elementos_nos"]["elementos_dominio"]
        
        # Extrair vetor de tempo (10 dias = 864000 s)
        vetor_tempo = dados["vetor_tempo"]
        
        print(f"Bloco 1 - Domínios ativos: {dominios_ativos}")
        print(f"Bloco 1 - Elementos ativos: {len(elementos_ativos_bloco1)} elementos")
        print(f"Bloco 1 - Nós de contorno: {len(nos_contorno)} nós")
        print(f"Vetor de tempo: {len(vetor_tempo)} pontos, até {vetor_tempo[-1]/86400:.1f} dias")
        
        # Criar submesh com elementos ativos do Bloco 1
        print(f"\nCriando submesh com elementos do Bloco 1...")
        elementos_ativos_array = np.array(elementos_ativos_bloco1, dtype=np.int32)
        
        from dolfinx import mesh as dolfinx_mesh
        submesh_result = dolfinx_mesh.create_submesh(mesh_global, mesh_global.topology.dim, elementos_ativos_array)
        domain = submesh_result[0]  # A submesh é o primeiro elemento retornado
        
        print(f"Submesh criada: {domain.topology.index_map(domain.topology.dim).size_global} células")
        print(f"Redução: {len(elementos_ativos_bloco1)/mesh_global.topology.index_map(mesh_global.topology.dim).size_global*100:.1f}% da malha global")
        
    except Exception as e:
        print(f"ERRO ao carregar JSON: {e}")
        return False
    
    # Definir espaços de funções na SUBMESH
    V_T = fem.functionspace(domain, ("Lagrange", 1))    # Temperatura
    V_teq = fem.functionspace(domain, ("Lagrange", 1))  # Tempo equivalente
    V_Q = fem.functionspace(domain, ("Lagrange", 1))    # Geração de calor
    
    print(f"\nEspaços de funções (submesh):")
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

    bc_T = fem.dirichletbc(T_bc_func, boundary_dofs)
    bcs = [bc_T]
    
    print(f"Condição de contorno Dirichlet aplicada: T = {T_bc} °C")
    
    # Definir funções de teste e trial para temperatura
    u = ufl.TrialFunction(V_T)
    v = ufl.TestFunction(V_T)
    
    # Medidas de integração na submesh
    dx = ufl.Measure("dx", domain=domain)
    
    print("\nDefinindo forma variacional para temperatura...")
    
    # A forma bilinear será definida dinamicamente no loop (pois dt varia)
    # Por enquanto, definimos apenas a parte de difusão que é constante
    a_difusao = k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    
    print("Forma variacional base definida")
    
    # Resolver sistema acoplado
    print("\n=== INICIANDO SIMULAÇÃO ACOPLADA ===")
    
    # Abrir arquivos XDMF para salvar evolução temporal (ETAPA 4)
    try:
        xdmf_temp = io.XDMFFile(MPI.COMM_WORLD, "etapa4_temperatura.xdmf", "w")
        xdmf_teq = io.XDMFFile(MPI.COMM_WORLD, "etapa4_teq.xdmf", "w")
        xdmf_Q = io.XDMFFile(MPI.COMM_WORLD, "etapa4_Q.xdmf", "w")
        
        # Salvar malha uma única vez
        xdmf_temp.write_mesh(domain)
        xdmf_teq.write_mesh(domain)
        xdmf_Q.write_mesh(domain)
        
        print("Arquivos XDMF preparados para salvamento temporal")
    except Exception as e:
        print(f"Erro ao preparar arquivos XDMF: {e}")
        return False
    
    # Executar simulação usando vetor de tempo do JSON (pula o primeiro ponto t=0)
    total_steps = len(vetor_tempo) - 1
    print(f"Total de passos de tempo: {total_steps}")
    
    for step in range(1, len(vetor_tempo)):
        t_anterior = vetor_tempo[step-1]
        t_atual = vetor_tempo[step]
        dt_atual = t_atual - t_anterior
        
        # Mostrar progresso apenas para alguns passos (evitar spam na tela)
        show_details = (step <= 5 or step % 20 == 0 or step == total_steps or t_atual % 86400 < dt_atual)
        
        if show_details:
            print(f"\nPasso {step}/{total_steps}:")
            print(f"  Tempo = {t_atual} s = {t_atual / 86400:.2f} dias")
            print(f"  dt = {dt_atual} s = {dt_atual / 3600:.2f} h")
        elif step % 10 == 0:
            print(f"  Passo {step}/{total_steps} - {t_atual / 86400:.1f} dias")
        
        # Trabalhando com submesh do Bloco 1 para toda a simulação
        
        # ETAPA 4.1: Atualizar tempo equivalente
        if show_details:
            print("  4.1 - Atualizando tempo equivalente...")
        teq_new = update_equivalent_time(teq_anterior.x.array, dt_atual, T_anterior.x.array, T_ref, A)
        teq_atual.x.array[:] = teq_new
        
        teq_min = np.min(teq_atual.x.array)
        teq_max = np.max(teq_atual.x.array)
        teq_mean = np.mean(teq_atual.x.array)
        if show_details:
            print(f"    teq_min = {teq_min:.2f} s")
            print(f"    teq_max = {teq_max:.2f} s")
            print(f"    teq_mean = {teq_mean:.2f} s")
        
        # ETAPA 4.2: Calcular geração de calor
        if show_details:
            print("  4.2 - Calculando geração de calor...")
        Q_new, alpha = calculate_heat_generation(teq_atual.x.array, alpha_max, tau, Q_total)
        Q_atual.x.array[:] = Q_new
        
        Q_min = np.min(Q_atual.x.array)
        Q_max_val = np.max(Q_atual.x.array)
        Q_mean = np.mean(Q_atual.x.array)
        alpha_mean = np.mean(alpha)
        if show_details:
            print(f"    Q_min = {Q_min:.2f} W/m³")
            print(f"    Q_max = {Q_max_val:.2f} W/m³")
            print(f"    Q_mean = {Q_mean:.2f} W/m³")
            print(f"    alpha_mean = {alpha_mean:.4f}")
        
        # ETAPA 4.3: Resolver equação de temperatura
        if show_details:
            print("  4.3 - Resolvendo equação de temperatura...")
        
        # Forma bilinear com dt dinâmico (termo transiente + difusão)
        a = (rho * ce / dt_atual) * u * v * ufl.dx + a_difusao
        
        # Forma linear com geração de calor atual
        L = (rho * ce / dt_atual) * T_anterior * v * ufl.dx
        L += Q_atual * v * ufl.dx  # Adicionar geração de calor
        
        # Usar LinearProblem que é mais robusto e funciona melhor
        from dolfinx.fem.petsc import LinearProblem
        
        try:
            # Criar e resolver problema linear
            problem = LinearProblem(a, L, bcs=bcs, 
                                   petsc_options={"ksp_type": "preonly", 
                                                 "pc_type": "lu"})
            
            # Resolver o problema - isso retorna uma nova função
            T_nova = problem.solve()
            
            # Copiar solução para T_atual
            T_atual.x.array[:] = T_nova.x.array[:]
            
            if show_details:
                print(f"    Solução obtida com sucesso")
            
            # Diagnósticos de temperatura
            T_min = np.min(T_atual.x.array)
            T_max = np.max(T_atual.x.array)
            T_mean = np.mean(T_atual.x.array)
            
            if show_details:
                print(f"    T_min = {T_min:.6f} °C")
                print(f"    T_max = {T_max:.6f} °C")
                print(f"    T_mean = {T_mean:.6f} °C")
            
            # Verificar se solução é finita
            if np.all(np.isfinite(T_atual.x.array)):
                if show_details:
                    print(f"    ✓ Solução é finita")
            else:
                print(f"    ✗ Solução contém valores não-finitos!")
                return False
            
            # Verificar se não divergiu para infinito
            if T_max < 1000.0:  # Limite razoável
                if show_details:
                    print(f"    ✓ Temperatura dentro de limites razoáveis")
            else:
                print(f"    ✗ Temperatura muito alta - possível divergência!")
                return False
            
            # Salvar resultados para este passo de tempo
            try:
                xdmf_temp.write_function(T_atual, t_atual)
                xdmf_teq.write_function(teq_atual, t_atual)
                xdmf_Q.write_function(Q_atual, t_atual)
                
                if show_details:
                    print(f"    ✓ Resultados salvos para t = {t_atual/86400:.2f} dias")
            except Exception as e_save:
                print(f"    ✗ Erro ao salvar passo {step}: {e_save}")
            
            # Atualizar campos anteriores
            T_anterior.x.array[:] = T_atual.x.array[:]
            teq_anterior.x.array[:] = teq_atual.x.array[:]
            
        except Exception as e:
            print(f"    ✗ Erro na solução: {e}")
            return False
    
    # Fechar arquivos XDMF
    print(f"\nFinalizando salvamento...")
    
    try:
        xdmf_temp.close()
        xdmf_teq.close()
        xdmf_Q.close()
        
        print("✓ Resultados salvos com sucesso")
        print(f"  - etapa4_temperatura.xdmf ({total_steps} passos de tempo)")
        print(f"  - etapa4_teq.xdmf ({total_steps} passos de tempo)")
        print(f"  - etapa4_Q.xdmf ({total_steps} passos de tempo)")
    except Exception as e:
        print(f"✗ Erro ao finalizar salvamento: {e}")
    
    print("\n=== ETAPA 4 CONCLUÍDA COM SUCESSO ===")
    return True

if __name__ == "__main__":
    success = main()
    if success:
        print("RESULTADO: SUCESSO")
        exit(0)
    else:
        print("RESULTADO: FALHA")
        exit(1) 