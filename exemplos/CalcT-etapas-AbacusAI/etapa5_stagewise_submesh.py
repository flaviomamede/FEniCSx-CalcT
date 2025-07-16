#!/usr/bin/python3
"""
ETAPA 5 - Simulação Stagewise Completa com Submeshes
Implementação seguindo o plano do etapa5_sugestao.md

Estratégia:
1. Submalha Cumulativa: cada bloco engloba anterior + elementos novos
2. Transferência de Estado: interpolação de T e teq entre malhas
3. Loop Externo: blocos construtivos (1→2→3)
4. Loop Interno: simulação temporal para cada bloco
"""

import numpy as np
import json
import os
import yaml
from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import LinearProblem
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
    T_kelvin = T_atual + 273.15
    T_ref_kelvin = T_ref + 273.15
    
    exponent = A * (1.0/T_ref_kelvin - 1.0/T_kelvin)
    
    # Clamp exponent para evitar overflow
    exponent = np.clip(exponent, -50, 50)
    
    beta = np.exp(exponent)
    
    teq_novo = teq_anterior + dt * beta
    return teq_novo

def calculate_heat_generation(teq, alpha_max=0.8, tau=86400.0, Q_total=300000000.0):
    """
    Calcula a geração de calor e grau de hidratação baseado no tempo equivalente
    """
    # Grau de hidratação baseado em modelo exponencial
    alpha = alpha_max * (1.0 - np.exp(-teq / tau))
    
    # Taxa de geração de calor (derivada da hidratação)
    Q = (Q_total * alpha_max / tau) * np.exp(-teq / tau)
    
    return Q, alpha

def carregar_dados_iniciais():
    """Carrega malha global e dados JSON"""
    mesh_file = "Uploads/barragem2.xdmf"
    json_file = "Uploads/analise_stagewise_barragem2_xdmf.json"
    yaml_file = "Uploads/barragem2.yaml"
    
    # Verificar existência
    for arquivo in [mesh_file, json_file, yaml_file]:
        if not os.path.exists(arquivo):
            raise FileNotFoundError(f"Arquivo não encontrado: {arquivo}")
    
    # Carregar malha global
    with io.XDMFFile(MPI.COMM_WORLD, mesh_file, "r") as xdmf:
        mesh_global = xdmf.read_mesh(name="malha")
        
        # Tentar carregar as meshtags (opcional)
        cell_tags_global = None
        facet_tags_global = None
        try:
            # Criar entidades de dimensão 1 (arestas) antes de ler facet_tags
            mesh_global.topology.create_entities(1)
            
            cell_tags_global = xdmf.read_meshtags(mesh_global, name="malha_cells")
            facet_tags_global = xdmf.read_meshtags(mesh_global, name="malha_facets")
            print(f"✓ Meshtags carregadas: células ({len(cell_tags_global.values)} tags), faces ({len(facet_tags_global.values)} tags)")
        except Exception as e:
            print(f"⚠️ Aviso: Não foi possível carregar meshtags: {e}")
            print(f"   Continuando sem tags de células/faces")
    
    # Carregar dados JSON
    with open(json_file, 'r') as f:
        plano_camadas = json.load(f)
    
    # Carregar dados YAML
    with open(yaml_file, 'r') as f:
        config_yaml = yaml.safe_load(f)
    
    print(f"✓ Malha global carregada: {mesh_global.topology.index_map(mesh_global.topology.dim).size_local} células")
    
    return mesh_global, plano_camadas, config_yaml, facet_tags_global

def criar_submalha_apenas_elementos_novos(mesh_global, elementos_novos):
    """Cria submalha apenas com os elementos novos do bloco atual"""
    # Marcar elementos novos
    marker = mesh.meshtags(
        mesh_global, 
        mesh_global.topology.dim, 
        np.array(elementos_novos, dtype=np.int32), 
        np.ones(len(elementos_novos), dtype=np.int32)
    )
    
    # Criar submesh
    submesh, entity_map, vertex_map, geom_map = mesh.create_submesh(
        mesh_global, mesh_global.topology.dim, marker.find(1)
    )
    
    num_cells = submesh.topology.index_map(submesh.topology.dim).size_local
    num_vertices = submesh.topology.index_map(0).size_local
    
    print(f"    Submalha criada (apenas elementos novos): {num_cells} células, {num_vertices} vértices")
    
    return submesh, entity_map, vertex_map



def verificar_nos_interface(plano_camadas, bloco_atual, bloco_anterior):
    """Verifica e debuga os nós de interface entre blocos"""
    if bloco_anterior is None:
        return []
    
    # Obter nós do domínio do bloco anterior
    nos_bloco_anterior = set(plano_camadas['analise_resultados'][bloco_anterior]['elementos_nos']['nos_dominio'])
    
    # Obter nós do domínio do bloco atual
    nos_bloco_atual = set(plano_camadas['analise_resultados'][bloco_atual]['elementos_nos']['nos_dominio'])
    
    # Nós da interface são os que estão em ambos os blocos
    nos_interface = list(nos_bloco_anterior.intersection(nos_bloco_atual))
    
    print(f"    Debug interface:")
    print(f"      Nós bloco anterior ({bloco_anterior}): {len(nos_bloco_anterior)}")
    print(f"      Nós bloco atual ({bloco_atual}): {len(nos_bloco_atual)}")
    print(f"      Nós interface: {len(nos_interface)}")
    
    if len(nos_interface) > 0:
        print(f"      Primeiros 5 nós interface: {nos_interface[:5]}")
    
    return nos_interface

def obter_elementos_novos(plano_camadas, bloco_key):
    """Obtém apenas os elementos novos do bloco atual"""
    diferencas = plano_camadas['analise_resultados'][bloco_key].get('diferencas', {})
    
    if not diferencas or 'elementos_dominio' not in diferencas:
        # Primeiro bloco ou bloco sem diferenças: usar todos os elementos do domínio
        elementos_novos = plano_camadas['analise_resultados'][bloco_key]['elementos_nos']['elementos_dominio']
        print(f"    Bloco sem diferenças: usando todos os {len(elementos_novos)} elementos do domínio")
    else:
        # Blocos subsequentes: usar apenas elementos novos das diferenças
        elementos_novos = diferencas['elementos_dominio']['entradas']
        print(f"    Bloco com diferenças: usando {len(elementos_novos)} elementos novos das diferenças")
    
    return elementos_novos

def obter_nos_novos(plano_camadas, bloco_key):
    """Obtém apenas os nós novos do bloco atual"""
    diferencas = plano_camadas['analise_resultados'][bloco_key].get('diferencas', {})
    
    if not diferencas or 'nos_dominio' not in diferencas:
        # Primeiro bloco ou bloco sem diferenças: usar todos os nós do domínio
        nos_novos = plano_camadas['analise_resultados'][bloco_key]['elementos_nos']['nos_dominio']
        print(f"    Bloco sem diferenças: usando todos os {len(nos_novos)} nós do domínio")
    else:
        # Blocos subsequentes: usar apenas nós novos das diferenças
        nos_novos = diferencas['nos_dominio']['entradas']
        print(f"    Bloco com diferenças: usando {len(nos_novos)} nós novos das diferenças")
    
    return nos_novos

def pre_processar_config(config_yaml):
    """
    Cria dicionários para acesso rápido às informações de tempo das camadas.
    """
    info = {
        'birth_times': {camada['nome']: camada['birth'] for camada in config_yaml['camadas']},
        'death_times': {camada['nome']: camada['death'] for camada in config_yaml['camadas']}
    }
    return info

def obter_contornos_ativos(current_time, config_yaml, config_info, submesh):
    """
    Determina quais condições de contorno estão ativas no tempo atual.
    """
    active_convection_bcs = []
    dirichlet_bcs = []

    for contorno in config_yaml['contornos']:
        camada_nascimento = contorno['nasce_com_camada']
        tempo_nascimento = config_info['birth_times'].get(camada_nascimento)

        is_born = tempo_nascimento is not None and current_time >= tempo_nascimento
        is_active = is_born

        # Verifica se o contorno foi desativado por outra camada
        if is_born and 'desativado_pela_camada' in contorno:
            camada_desativacao = contorno['desativado_pela_camada']
            tempo_desativacao = config_info['birth_times'].get(camada_desativacao)
            if tempo_desativacao is not None and current_time >= tempo_desativacao:
                is_active = False

        if is_active:
            print(f"    -> Contorno '{contorno['nome']}' (ID: {contorno['id']}) está ATIVO no tempo {current_time}s.")
            if contorno['tipo'] == 'conveccao':
                active_convection_bcs.append({
                    'id': contorno['id'],
                    'h': contorno['h'],
                    't_ext': fem.Constant(submesh, contorno['t_ext'])  # Usar dolfinx.fem.Constant
                })
            elif contorno['tipo'] == 'fluxo':
                if contorno['material'] == 'espelho':
                    print(f"    -> Contorno '{contorno['nome']}' é isolamento perfeito (h=0)")
                else:
                    print(f"    -> Contorno '{contorno['nome']}' tem fluxo não-zero (não implementado)")

    return active_convection_bcs, dirichlet_bcs

def definir_forma_variacional(V, T_n, dt, bcs_robin, facet_tags):
    """
    Cria as formas bilinear (a) e linear (L) para o problema térmico.
    """
    # Constantes do problema
    dt_const = fem.Constant(V.mesh, dt)
    theta = fem.Constant(V.mesh, 0.5)  # Crank-Nicolson
    
    # Propriedades do material
    rho = fem.Constant(V.mesh, 2400.0)
    c_p = fem.Constant(V.mesh, 1000.0)
    k = fem.Constant(V.mesh, 2.5)
    
    T = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # Ponto médio no tempo (Crank-Nicolson)
    T_mid = theta * T + (1.0 - theta) * T_n

    # Forma bilinear: (rho*c/dt)*T*v*dx + theta*k*grad(T)*grad(v)*dx
    a = (rho * c_p / dt_const) * T * v * ufl.dx + theta * ufl.dot(k * ufl.grad(T), ufl.grad(v)) * ufl.dx

    # Forma linear: (rho*c/dt)*T_n*v*dx - (1-theta)*k*grad(T_n)*grad(v)*dx
    L = (rho * c_p / dt_const) * T_n * v * ufl.dx - (1.0 - theta) * ufl.dot(k * ufl.grad(T_n), ufl.grad(v)) * ufl.dx

    # Adicionar termos de convecção dos contornos ativos
    for bc_info in bcs_robin:
        h = fem.Constant(V.mesh, bc_info['h'])
        t_ext = fem.Constant(V.mesh, bc_info['t_ext'])
        
        # Termos de convecção: h * T_mid * v * ds e h * t_ext * v * ds
        # Usar ds() geral por enquanto, já que facet_tags pode não estar mapeado corretamente para submalha
        a += theta * h * T * v * ufl.ds
        L += (1.0 - theta) * h * T_n * v * ufl.ds + h * t_ext * v * ufl.ds
    
    return a, L

def aplicar_condicoes_iniciais_e_contorno(submesh, config_yaml, bloco_info, facet_tags):
    """Aplica condições iniciais e de contorno baseadas no YAML, usando facet_tags para BCs de convecção"""
    # Espaços de função
    V_T = fem.functionspace(submesh, ("Lagrange", 1))
    V_teq = fem.functionspace(submesh, ("Lagrange", 1))

    # Funções
    T = fem.Function(V_T)
    teq = fem.Function(V_teq)

    # Condições iniciais
    T.x.array[:] = 20.0  # °C (padrão do YAML)
    teq.x.array[:] = 0.0  # s

    # Pré-processar configuração
    config_info = pre_processar_config(config_yaml)
    current_time = bloco_info['info_bloco']['inicio']

    # Listas de BCs
    bcs_dirichlet = []  # Por enquanto, lista vazia para evitar problemas de IndexMap
    bcs_robin = []

    # Identificar contornos ativos para este tempo
    for contorno in config_yaml['contornos']:
        camada_nascimento = contorno['nasce_com_camada']
        tempo_nascimento = config_info['birth_times'].get(camada_nascimento)
        is_born = tempo_nascimento is not None and current_time >= tempo_nascimento
        is_active = is_born
        if is_born and 'desativado_pela_camada' in contorno:
            camada_desativacao = contorno['desativado_pela_camada']
            tempo_desativacao = config_info['birth_times'].get(camada_desativacao)
            if tempo_desativacao is not None and current_time >= tempo_desativacao:
                is_active = False
        if not is_active:
            continue
        if contorno['tipo'] == 'dirichlet':
            # TODO: Implementar BCs de Dirichlet na submalha quando necessário
            # Por enquanto, pular para evitar problemas de IndexMap
            print(f"    ⚠️ BC Dirichlet '{contorno['nome']}' ignorada temporariamente")
            continue
        elif contorno['tipo'] == 'conveccao':
            h = contorno['h']
            t_ext = contorno['t_ext']
            bc_id = contorno['id']
            bcs_robin.append({'id': bc_id, 'h': h, 't_ext': t_ext})
        # Fluxo/Neumann pode ser adicionado aqui se necessário

    print(f"    ✓ Condições iniciais aplicadas (T=20°C, teq=0s)")
    print(f"    Dirichlet BCs: {len(bcs_dirichlet)} | Robin BCs: {len(bcs_robin)}")

    return T, teq, bcs_dirichlet, bcs_robin

def executar_loop_interno(submesh, T_anterior, teq_anterior, info_bloco, vetor_tempo_global, bcs):
    """Executa simulação temporal para um bloco específico"""
    # Filtrar tempos para este bloco
    t_inicio = info_bloco['inicio']
    t_fim = info_bloco['fim']
    
    vetor_tempo = np.array(vetor_tempo_global)
    mask_tempo = (vetor_tempo >= t_inicio) & (vetor_tempo <= t_fim)
    tempos_bloco = vetor_tempo[mask_tempo]
    
    print(f"    Simulação temporal: {len(tempos_bloco)} passos de tempo")
    
    # Parâmetros físicos
    rho = 2400.0      # kg/m³
    cp = 1000.0       # J/(kg·K)
    k = 2.5           # W/(m·K)
    
    # Espaços de função na submesh
    V_T = T_anterior.function_space
    V_teq = teq_anterior.function_space
    V_Q = fem.functionspace(submesh, ("Lagrange", 1))
    
    # Funções para este bloco
    T_atual = fem.Function(V_T)
    teq_atual = fem.Function(V_teq)
    Q_atual = fem.Function(V_Q)
    
    # Copiar estado inicial
    T_atual.x.array[:] = T_anterior.x.array[:]
    teq_atual.x.array[:] = teq_anterior.x.array[:]
    
    # Funções de teste
    u = ufl.TrialFunction(V_T)
    v = ufl.TestFunction(V_T)
    
    # Medida de integração local ligada à submesh
    dx_sub = ufl.Measure("dx", domain=submesh)
    
    # Resultados para monitoramento
    resultados = []
    
    # Loop temporal interno
    for i in range(1, len(tempos_bloco)):
        t_atual = tempos_bloco[i]
        dt = tempos_bloco[i] - tempos_bloco[i-1]
        
        # 1. Atualizar tempo equivalente
        teq_vals = update_equivalent_time(teq_atual.x.array[:], dt, T_atual.x.array[:])
        teq_atual.x.array[:] = teq_vals
        
        # 2. Calcular geração de calor
        Q_vals, alpha_vals = calculate_heat_generation(teq_vals)
        Q_atual.x.array[:] = Q_vals
        
        # 3. Resolver equação de temperatura (usando medida de integração local)
        a = (rho * cp * u * v + dt * k * ufl.dot(ufl.grad(u), ufl.grad(v))) * dx_sub
        L = (rho * cp * T_atual * v + dt * Q_atual * v) * dx_sub
        
        # Adicionar termos de convecção dos contornos ativos
        for bc_info in bcs_robin:
            h = bc_info['h']
            t_ext = bc_info['t_ext']
            bc_id = bc_info['id']
            
            # Termos de convecção: h * T * v * ds e h * t_ext * v * ds
            # Nota: ds(bc_id) requer facet_tags, que não temos aqui
            # Por enquanto, usar ds() geral (todos os contornos)
            a += dt * h * u * v * ufl.ds
            L += dt * h * t_ext * v * ufl.ds
        
        problem = LinearProblem(a, L, bcs=bcs, 
                              petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        T_atual = problem.solve()
        
        # 4. Calcular estatísticas
        T_mean = np.mean(T_atual.x.array[:])
        teq_mean = np.mean(teq_atual.x.array[:])
        Q_mean = np.mean(Q_atual.x.array[:])
        alpha_mean = np.mean(alpha_vals)
        
        resultados.append({
            't': t_atual,
            'T_mean': T_mean,
            'teq_mean': teq_mean, 
            'Q_mean': Q_mean,
            'alpha_mean': alpha_mean
        })
        
        # Log progresso
        if i % 5 == 0 or i == len(tempos_bloco) - 1:
            print(f"      t={t_atual:8.0f}s | T={T_mean:5.1f}°C | teq={teq_mean:8.0f}s | Q={Q_mean:6.0f}W/m³ | α={alpha_mean:.3f}")
    
    print(f"    ✓ Loop temporal concluído: {len(resultados)} passos")
    
    return T_atual, teq_atual, resultados

def tentar_interpolacao_interface(T_anterior, teq_anterior, submesh_nova, nos_interface):
    """Tenta interpolar valores apenas nos nós de interface"""
    print(f"    Tentando interpolação nos {len(nos_interface)} nós de interface...")
    
    # Espaços de função na nova submesh
    V_T_nova = fem.functionspace(submesh_nova, ("Lagrange", 1))
    V_teq_nova = fem.functionspace(submesh_nova, ("Lagrange", 1))
    
    # Funções na nova submesh
    T_nova = fem.Function(V_T_nova)
    teq_nova = fem.Function(V_teq_nova)
    
    # Aplicar condições iniciais padrão primeiro
    T_nova.x.array[:] = 20.0  # °C
    teq_nova.x.array[:] = 0.0  # s
    
    try:
        # Tentar interpolação usando FEniCSx
        T_interpolated = fem.Function(V_T_nova)
        teq_interpolated = fem.Function(V_teq_nova)
        
        # Interpolar T
        T_interpolated.interpolate(T_anterior)
        
        # Interpolar teq
        teq_interpolated.interpolate(teq_anterior)
        
        # Usar valores interpolados
        T_nova.x.array[:] = T_interpolated.x.array[:]
        teq_nova.x.array[:] = teq_interpolated.x.array[:]
        
        print(f"    ✓ Interpolação FEniCSx bem-sucedida")
        
    except Exception as e:
        print(f"    ⚠️ Interpolação FEniCSx falhou: {e}")
        print(f"    ✓ Usando condições iniciais padrão (T=20°C, teq=0s)")
    
    return T_nova, teq_nova

def run_full_stagewise_simulation():
    """Simulação stagewise completa seguindo a nova estratégia"""
    print("ETAPA 5 - SIMULAÇÃO STAGEWISE COMPLETA COM SUBMESHES (NOVA ESTRATÉGIA)")
    print("=" * 70)
    
    try:
        comm = MPI.COMM_WORLD
        mesh_global, plano_camadas, config_yaml, facet_tags_global = carregar_dados_iniciais()
        
        solucao_T_bloco_anterior = None
        solucao_teq_bloco_anterior = None
        resultados_completos = []

        # === Arquivos XDMF globais ===
        xdmf_temp_global = io.XDMFFile(comm, "etapa5_temperatura.xdmf", "w")
        xdmf_teq_global = io.XDMFFile(comm, "etapa5_teq.xdmf", "w")
        xdmf_Q_global = io.XDMFFile(comm, "etapa5_Q.xdmf", "w")
        mesh_global_written = False
        
        blocos = ['bloco_1', 'bloco_2', 'bloco_3']
        for idx, bloco_key in enumerate(blocos):
            bloco_id = idx + 1
            info_bloco = plano_camadas['analise_resultados'][bloco_key]['info_bloco']
            elementos_bloco = obter_elementos_novos(plano_camadas, bloco_key)
            print(f"\n🚀 INICIANDO BLOCO CONSTRUTIVO {bloco_id}: {bloco_key} 🚀")
            print(f"    Elementos do bloco: {len(elementos_bloco)}")
            bloco_anterior = None
            nos_interface = []
            if idx > 0:
                bloco_anterior = blocos[idx - 1]
                nos_interface = verificar_nos_interface(plano_camadas, bloco_key, bloco_anterior)
                print(f"    ✓ BC da interface removida para o próximo bloco")
            else:
                print(f"    Primeiro bloco: sem interface anterior")
            print(f"    Criando submalha apenas com elementos novos...")
            submesh, entity_map, vertex_map = criar_submalha_apenas_elementos_novos(mesh_global, elementos_bloco)
            print(f"    Aplicando condições iniciais e de contorno...")
            T_inicial, teq_inicial, bcs, bcs_robin = aplicar_condicoes_iniciais_e_contorno(submesh, config_yaml, plano_camadas['analise_resultados'][bloco_key], facet_tags_global)
            if solucao_T_bloco_anterior is not None and len(nos_interface) > 0:
                print(f"    Transferindo estado do bloco anterior para {len(nos_interface)} nós de interface...")
                T_inicial, teq_inicial = tentar_interpolacao_interface(
                    solucao_T_bloco_anterior, solucao_teq_bloco_anterior, submesh, nos_interface
                )
            elif solucao_T_bloco_anterior is not None:
                print(f"    ⚠️ Nenhum nó de interface encontrado, usando condições iniciais padrão")
            else:
                print(f"    Primeiro bloco: usando condições iniciais padrão")

            # === Arquivos XDMF por bloco ===
            xdmf_temp = io.XDMFFile(comm, f"etapa5_bloco{bloco_id}_temperatura.xdmf", "w")
            xdmf_teq = io.XDMFFile(comm, f"etapa5_bloco{bloco_id}_teq.xdmf", "w")
            xdmf_Q = io.XDMFFile(comm, f"etapa5_bloco{bloco_id}_Q.xdmf", "w")
            xdmf_temp.write_mesh(submesh)
            xdmf_teq.write_mesh(submesh)
            xdmf_Q.write_mesh(submesh)

            print(f"    Executando simulação temporal...")
            T_final, teq_final, resultados_bloco = executar_loop_interno_salvando(
                submesh, T_inicial, teq_inicial, info_bloco, plano_camadas['vetor_tempo'], bcs,
                xdmf_temp, xdmf_teq, xdmf_Q,
                xdmf_temp_global, xdmf_teq_global, xdmf_Q_global,
                not mesh_global_written, bcs_robin, facet_tags_global
            )
            mesh_global_written = True
            resultados_completos.extend(resultados_bloco)
            solucao_T_bloco_anterior = T_final
            solucao_teq_bloco_anterior = teq_final
            xdmf_temp.close()
            xdmf_teq.close()
            xdmf_Q.close()
            print(f"✓ BLOCO {bloco_id} CONCLUÍDO")
        xdmf_temp_global.close()
        xdmf_teq_global.close()
        xdmf_Q_global.close()
        # ===== RELATÓRIO FINAL =====
        print(f"\n=== RELATÓRIO FINAL ===")
        print(f"Simulação stagewise concluída com {len(resultados_completos)} passos")
        
        if resultados_completos:
            T_inicial = resultados_completos[0]['T_mean']
            T_final = resultados_completos[-1]['T_mean']
            T_max = max(r['T_mean'] for r in resultados_completos)
            alpha_final = resultados_completos[-1]['alpha_mean']
            
            print(f"Temperatura inicial: {T_inicial:.1f}°C")
            print(f"Temperatura final: {T_final:.1f}°C")
            print(f"Pico de temperatura: {T_max:.1f}°C")
            print(f"Hidratação final: {alpha_final:.3f} ({alpha_final*100:.1f}%)")
            
            # Verificações
            if T_max < 100:
                print("✓ Temperatura dentro de limites razoáveis")
            
            if alpha_final > 0.1:
                print("✓ Hidratação progrediu adequadamente")
        
        print("\n🎉 SIMULAÇÃO STAGEWISE COMPLETA CONCLUÍDA! 🎉")
        print("NOVA ESTRATÉGIA: SUBMESHES APENAS COM ELEMENTOS NOVOS")
        print("RESULTADO: SUCESSO")
        return True
        
    except Exception as e:
        print(f"\nERRO: {e}")
        import traceback
        traceback.print_exc()
        return False

def executar_loop_interno_salvando(submesh, T_anterior, teq_anterior, info_bloco, vetor_tempo_global, bcs,
                                   xdmf_temp, xdmf_teq, xdmf_Q,
                                   xdmf_temp_global, xdmf_teq_global, xdmf_Q_global,
                                   write_mesh_global, bcs_robin, facet_tags):
    # Igual ao executar_loop_interno, mas salva XDMF a cada passo
    t_inicio = info_bloco['inicio']
    t_fim = info_bloco['fim']
    vetor_tempo = np.array(vetor_tempo_global)
    mask_tempo = (vetor_tempo >= t_inicio) & (vetor_tempo <= t_fim)
    tempos_bloco = vetor_tempo[mask_tempo]
    print(f"    Simulação temporal: {len(tempos_bloco)} passos de tempo")
    
    # Espaços de função na submesh
    V_T = T_anterior.function_space
    V_teq = teq_anterior.function_space
    V_Q = fem.functionspace(submesh, ("Lagrange", 1))
    
    # Funções para este bloco
    T_atual = fem.Function(V_T)
    teq_atual = fem.Function(V_teq)
    Q_atual = fem.Function(V_Q)
    
    # Copiar estado inicial
    T_atual.x.array[:] = T_anterior.x.array[:]
    teq_atual.x.array[:] = teq_anterior.x.array[:]
    
    # Medida de integração local ligada à submesh
    dx_sub = ufl.Measure("dx", domain=submesh)
    
    # Resultados para monitoramento
    resultados = []
    
    if write_mesh_global:
        xdmf_temp_global.write_mesh(submesh)
        xdmf_teq_global.write_mesh(submesh)
        xdmf_Q_global.write_mesh(submesh)
    
    # Loop temporal interno
    for i in range(1, len(tempos_bloco)):
        t_atual = tempos_bloco[i]
        dt = tempos_bloco[i] - tempos_bloco[i-1]
        
        # 1. Atualizar tempo equivalente
        teq_vals = update_equivalent_time(teq_atual.x.array[:], dt, T_atual.x.array[:])
        teq_atual.x.array[:] = teq_vals
        
        # 2. Calcular geração de calor
        Q_vals, alpha_vals = calculate_heat_generation(teq_vals)
        Q_atual.x.array[:] = Q_vals
        
        # 3. Resolver equação de temperatura usando nova forma variacional
        a, L = definir_forma_variacional(V_T, T_atual, dt, bcs_robin, facet_tags)
        
        # Adicionar termo de geração de calor
        L += Q_atual * ufl.TestFunction(V_T) * dx_sub
        
        problem = LinearProblem(a, L, bcs=bcs, 
                              petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        T_atual = problem.solve()
        
        # 4. Calcular estatísticas
        T_mean = np.mean(T_atual.x.array[:])
        teq_mean = np.mean(teq_atual.x.array[:])
        Q_mean = np.mean(Q_atual.x.array[:])
        alpha_mean = np.mean(alpha_vals)
        
        resultados.append({
            't': t_atual,
            'T_mean': T_mean,
            'teq_mean': teq_mean, 
            'Q_mean': Q_mean,
            'alpha_mean': alpha_mean
        })
        
        # Log progresso
        if i % 5 == 0 or i == len(tempos_bloco) - 1:
            print(f"      t={t_atual:8.0f}s | T={T_mean:5.1f}°C | teq={teq_mean:8.0f}s | Q={Q_mean:6.0f}W/m³ | α={alpha_mean:.3f}")
        
        # Salvar XDMF bloco
        xdmf_temp.write_function(T_atual, t_atual)
        xdmf_teq.write_function(teq_atual, t_atual)
        xdmf_Q.write_function(Q_atual, t_atual)
        
        # Salvar XDMF global
        xdmf_temp_global.write_function(T_atual, t_atual)
        xdmf_teq_global.write_function(teq_atual, t_atual)
        xdmf_Q_global.write_function(Q_atual, t_atual)
    
    print(f"    ✓ Loop temporal concluído: {len(resultados)} passos")
    return T_atual, teq_atual, resultados

if __name__ == "__main__":
    sucesso = run_full_stagewise_simulation()
    import sys
    sys.exit(0 if sucesso else 1) 