#!/usr/bin/python3
"""
TESTE ETAPA 1 - Correção da Medida de Integração
Testa se o uso de dx_sub resolve o segmentation fault
Executa apenas o Bloco 1 para validação isolada
"""

import numpy as np
import json
import os
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import petsc4py
petsc4py.init()
from petsc4py import PETSc

def test_bloco1_medida_integracao():
    """Testa apenas o Bloco 1 com medida de integração corrigida"""
    print("TESTE ETAPA 1 - CORREÇÃO DA MEDIDA DE INTEGRAÇÃO")
    print("=" * 50)
    
    try:
        # ===== CARREGAR DADOS =====
        mesh_file = "Uploads/barragem2.xdmf"
        json_file = "Uploads/analise_stagewise_barragem2_xdmf.json"
        
        # Verificar arquivos
        for arquivo in [mesh_file, json_file]:
            if not os.path.exists(arquivo):
                raise FileNotFoundError(f"Arquivo não encontrado: {arquivo}")
        
        # Carregar malha global
        with io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
            mesh_global = xdmf.read_mesh(name="malha")
            cell_tags_global = xdmf.read_meshtags(mesh_global, name="malha_cells")
        
        # Carregar dados JSON
        with open(json_file, 'r') as f:
            dados_json = json.load(f)
        
        print(f"✓ Dados carregados: {mesh_global.topology.index_map(mesh_global.topology.dim).size_local} células")
        
        # ===== CRIAR SUBMESH DO BLOCO 1 =====
        bloco1_elementos = dados_json['analise_resultados']['bloco_1']['elementos_nos']['elementos_dominio']
        
        print(f"✓ Bloco 1: {len(bloco1_elementos)} elementos")
        
        # Marcar elementos do Bloco 1
        marker = mesh.meshtags(
            mesh_global, 
            mesh_global.topology.dim, 
            np.array(bloco1_elementos, dtype=np.int32), 
            np.ones(len(bloco1_elementos), dtype=np.int32)
        )
        
        # Criar submesh
        submesh, entity_map, vertex_map, geom_map = mesh.create_submesh(
            mesh_global, mesh_global.topology.dim, marker.find(1)
        )
        
        num_cells = submesh.topology.index_map(submesh.topology.dim).size_local
        num_vertices = submesh.topology.index_map(0).size_local
        
        print(f"✓ Submesh criada: {num_cells} células, {num_vertices} vértices")
        
        # ===== DEFINIR ESPAÇOS E FUNÇÕES =====
        V_T = fem.functionspace(submesh, ("Lagrange", 1))
        V_teq = fem.functionspace(submesh, ("Lagrange", 1))
        V_Q = fem.functionspace(submesh, ("Lagrange", 1))
        
        print(f"✓ Espaços definidos: {V_T.dofmap.index_map.size_local} DOFs")
        
        # Funções
        T_atual = fem.Function(V_T)
        T_anterior = fem.Function(V_T)
        teq_atual = fem.Function(V_teq)
        teq_anterior = fem.Function(V_teq)
        Q_atual = fem.Function(V_Q)
        
        # Condições iniciais
        T_atual.x.array[:] = 20.0    # °C
        T_anterior.x.array[:] = 20.0
        teq_atual.x.array[:] = 0.0   # s
        teq_anterior.x.array[:] = 0.0
        Q_atual.x.array[:] = 0.0     # W/m³
        
        print("✓ Funções e condições iniciais definidas")
        
        # ===== DEFINIR FORMA VARIACIONAL COM MEDIDA CORRETA =====
        u = ufl.TrialFunction(V_T)
        v = ufl.TestFunction(V_T)
        
        # CORREÇÃO CRÍTICA: Medida de integração local
        dx_sub = ufl.Measure("dx", domain=submesh)
        
        print("✓ Medida de integração local criada: dx_sub")
        
        # Parâmetros físicos
        rho = 2400.0      # kg/m³
        cp = 1000.0       # J/(kg·K)
        k = 2.5           # W/(m·K)
        dt = 3600.0       # s (1 hora para teste)
        
        # Condições de contorno (isolamento)
        bcs = []
        
        # ===== TESTAR MONTAGEM E RESOLUÇÃO =====
        print("Testando montagem da forma variacional...")
        
        # Forma bilinear e linear com medida local
        a = (rho * cp * u * v + dt * k * ufl.dot(ufl.grad(u), ufl.grad(v))) * dx_sub
        L = (rho * cp * T_anterior * v + dt * Q_atual * v) * dx_sub
        
        print("✓ Formas variacionais definidas com dx_sub")
        
        # Tentar resolver sistema
        print("Testando solução do sistema linear...")
        
        problem = LinearProblem(a, L, bcs=bcs, 
                              petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        T_atual = problem.solve()
        
        print("✓ Sistema linear resolvido com sucesso!")
        
        # Verificar resultado
        T_mean = np.mean(T_atual.x.array[:])
        print(f"✓ Temperatura média: {T_mean:.2f}°C")
        
        # ===== TESTAR ALGUNS PASSOS TEMPORAIS =====
        print("Testando múltiplos passos temporais...")
        
        for passo in range(1, 4):  # 3 passos para teste
            print(f"  Passo {passo}:")
            
            # Atualizar tempo equivalente (simplificado)
            teq_atual.x.array[:] = teq_anterior.x.array[:] + dt
            
            # Gerar calor (simplificado)
            Q_atual.x.array[:] = 1000.0  # W/m³ constante para teste
            
            # Resolver
            L = (rho * cp * T_anterior * v + dt * Q_atual * v) * dx_sub
            problem = LinearProblem(a, L, bcs=bcs, 
                                  petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
            T_atual = problem.solve()
            
            # Atualizar para próximo passo
            T_anterior.x.array[:] = T_atual.x.array[:]
            teq_anterior.x.array[:] = teq_atual.x.array[:]
            
            T_mean = np.mean(T_atual.x.array[:])
            print(f"    T_média = {T_mean:.2f}°C")
        
        print("\n🎉 TESTE ETAPA 1: SUCESSO!")
        print("✓ Medida de integração dx_sub funcionando corretamente")
        print("✓ Sem segmentation fault detectado")
        print("✓ Submesh do Bloco 1 operacional")
        
        return True
        
    except Exception as e:
        print(f"\n❌ ERRO: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    sucesso = test_bloco1_medida_integracao()
    import sys
    
    if sucesso:
        print("\n✅ ETAPA 1 CONCLUÍDA - Medida de integração corrigida!")
        print("📈 Próximo passo: testar interpolação entre blocos")
        sys.exit(0)
    else:
        print("\n❌ ETAPA 1 FALHADA - Verificar implementação")
        sys.exit(1) 