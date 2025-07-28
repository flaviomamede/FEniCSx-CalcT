#!/usr/bin/env python3
"""
BARRAGEM FEniCSx - VERSÃO FINAL REFATORADA (R2)

Este script implementa uma análise de transferência de calor em uma barragem
construída em etapas (stagewise), utilizando uma abordagem robusta e modular.

NOVIDADES DESTA VERSÃO (R2):
- **CORREÇÃO DE SINTAXE:** Corrigida a chamada para a criação de FunctionSpace
  para `dolfinx.fem.functionspace` (minúsculo), compatível com as versões
  recentes do FEniCSx.
- **Atribuição Correta de Múltiplos Materiais:** Utiliza funções (DG-0) para
  atribuir propriedades distintas para cada grupo de elementos.
"""

import os
import sys
import json
import yaml
import datetime
import numpy as np
import ufl
from dolfinx import mesh, io
from dolfinx.fem import (Function, Constant, locate_dofs_topological,
                         functionspace, form, dirichletbc)
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import petsc4py
from pathlib import Path

petsc4py.init()
from petsc4py import PETSc

class TeeOutput:
    """Classe para redirecionar a saída para o terminal e um arquivo de log."""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log_file = open(filename, 'w', encoding='utf-8')
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = f"# LOG DE EXECUÇÃO - SIMULAÇÃO DE BARRAGEM\n# Data: {timestamp}\n\n"
        self.log_file.write(header)
        self.log_file.flush()

    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
        self.log_file.flush()

    def flush(self):
        self.terminal.flush()
        self.log_file.flush()

    def close(self):
        self.log_file.close()
        sys.stdout = self.terminal

class SimulacaoBarragemR2:
    def __init__(self, config_file, json_file, log_file="log_simulacao.md"):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.config_file = config_file
        self.json_file = json_file
        self.log_file = log_file
        self._load_config_and_analysis()

    def _load_config_and_analysis(self):
        """Carrega o YAML de configuração e o JSON de análise pré-processado."""
        with open(self.config_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        with open(self.json_file, 'r', encoding='utf-8') as f:
            self.analysis = json.load(f)

        self.theta = self.config['general']['theta']
        
        # Define o diretório de saída dentro da pasta do caso
        case_dir = Path(self.config_file).parent
        self.output_dir = case_dir / self.config['general']['output_dir']
        
        mesh_file_rel = self.config['general']['mesh_file']
        self.mesh_file_path = case_dir / mesh_file_rel

        self.time_vector = self.analysis['vetor_tempo']
        self.time_blocks = self.analysis['blocos_tempo']
        self.simulation_plan = self.analysis['analise_resultados']
        self.mappings = self.analysis['info_geral']['mapeamentos']

        # Inicialização das temperaturas iniciais por camada/material
        self.temp_inicial = 20.0
        self.temp_iniciais_por_camada = {
            self.mappings['physical_surfaces'][camada_material['nome']]: 
            camada_material.get('temperatura_inicial', 20.0)
            for camada_material in self.config['camadas_material']
        }
       
    def run(self):
        try:
            self._setup()
            self._run_simulation_loop()
            self._finalize()
        except Exception as e:
            if self.rank == 0:
                print(f"\n❌ ERRO FATAL DURANTE A EXECUÇÃO: {e}")
                import traceback
                traceback.print_exc()
        finally:
            if hasattr(self, 'logger') and self.logger:
                self.logger.close()

    def _setup(self):
        if self.rank == 0:
            print("\n--- FASE DE CONFIGURAÇÃO ---")
        self._load_mesh()
        self._setup_materials_data()
        self._setup_function_spaces()
        self._assign_materials_to_mesh()
        os.makedirs(self.output_dir, exist_ok=True)

    def _load_mesh(self):
        if self.rank == 0:
            print(f"   ➡️  Carregando malha de '{self.mesh_file_path}'...")
        
        with io.XDMFFile(self.comm, str(self.mesh_file_path), "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            self.mesh.topology.create_entities(1)
            self.mesh.topology.create_connectivity(1, 2)
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        
        if self.rank == 0:
            print(f"   ✅ Malha '{self.mesh.name}' carregada.")

    def _setup_materials_data(self):
        """Carrega os dados dos materiais do YAML para um dicionário em memória."""
        if self.rank == 0:
            print("   ➡️  Carregando dados dos materiais do YAML...")
        self.materials = {}
        for mat in self.config['materiais']:
            hgen = mat.get('hgen', {})
            self.materials[mat['nome']] = {
                'densidade': mat['densidade'],
                'condutividade': mat['condutividade_termica'],
                'calor_especifico': mat['calor_especifico'],
                'gera_calor': hgen.get('gera_calor', False),
                'dTadinfty': hgen.get('par_gera_calor', {}).get('dTadinfty', 30.0),
                'a_dias': hgen.get('par_gera_calor', {}).get('a_dias', 1.5),
                'expoente': hgen.get('par_gera_calor', {}).get('expoente', 2.0),
                'EaR': hgen.get('EaR', 4000.0),
                'Tref': hgen.get('Tref', 20.0)
            }
        self.has_exothermic = any(m['gera_calor'] for m in self.materials.values())

    def _setup_function_spaces(self):
        """Define os espaços de função necessários para o problema."""
        if self.rank == 0:
            print("   ➡️  Definindo espaços de função...")
        
        self.V = functionspace(self.mesh, ("CG", 1))
        self.V_prop = functionspace(self.mesh, ("DG", 0))
        
        # Funções para armazenar as propriedades que variam espacialmente
        self.k = Function(self.V_prop, name="Condutividade")
        self.rho = Function(self.V_prop, name="Densidade")
        self.cp = Function(self.V_prop, name="CalorEspecifico")

        if self.has_exothermic:
            self.Tp = Function(self.V, name="Temperatura")
            self.Tp_n = Function(self.V)
            self.u_Tp = ufl.TrialFunction(self.V)
            self.v_Tp = ufl.TestFunction(self.V)
            self.teq = Function(self.V, name="TempoEquivalente")
            self.teq_n = Function(self.V)
            self.Q_heat = Function(self.V, name="GeracaoCalor")
        else:
            self.T = Function(self.V, name="Temperatura")
            self.T_n = Function(self.V)
            self.u = ufl.TrialFunction(self.V)
            self.v = ufl.TestFunction(self.V)
        
        if self.rank == 0:
            print("   ✅ Espaços de função definidos.")

    def _assign_materials_to_mesh(self):
        """Atribui as propriedades dos materiais aos elementos da malha (células)."""
        if self.rank == 0:
            print("   ➡️  Atribuindo propriedades dos materiais aos elementos da malha...")
        
        for camada_mat_info in self.config.get('camadas_material', []):
            nome_material = camada_mat_info.get('material')
            nome_pg_surf = camada_mat_info.get('nome')

            if nome_material in self.materials and nome_pg_surf in self.mappings['physical_surfaces']:
                pg_id = self.mappings['physical_surfaces'][nome_pg_surf]
                cells = self.cell_tags.find(pg_id)
                props = self.materials[nome_material]
                if cells.size > 0:
                    self.k.x.array[cells] = props['condutividade']
                    self.rho.x.array[cells] = props['densidade']
                    self.cp.x.array[cells] = props['calor_especifico']
                    if self.rank == 0:
                        print(f"      - Material '{nome_material}' atribuído a {len(cells)} elementos no PG {pg_id}.")

        # Cria arquivo de propriedades dentro da pasta do caso
        case_dir = Path(self.config_file).parent
        caso_name = case_dir.name
        prop_file = case_dir / f"{caso_name}-mat.xdmf"
        with io.XDMFFile(self.comm, prop_file, "w") as xdmf:
            xdmf.write_mesh(self.mesh)
            xdmf.write_function(self.k)
            xdmf.write_function(self.rho)
            xdmf.write_function(self.cp)
        if self.rank == 0:
            print(f"   ✅ Propriedades atribuídas. Verificação em: '{prop_file}'")

    def _run_simulation_loop(self):
        if self.rank == 0: print("\n--- FASE DE SIMULAÇÃO ---")
        total_steps = 0
        blocos = list(self.simulation_plan.items())
        for block_idx, (bloco_nome, bloco_info) in enumerate(blocos):
            block_time_points = [t for t in self.time_vector if bloco_info['info_bloco']['inicio'] <= t <= bloco_info['info_bloco']['fim']]
            if self.rank == 0:
                print("\n" + "="*80)
                print(f"📦 PROCESSANDO {bloco_nome.upper().replace('_', ' ')}")
                print(f"   - Período: {bloco_info['info_bloco']['inicio']/3600:.1f}h a {bloco_info['info_bloco']['fim']/3600:.1f}h")
                print(f"   - Pontos temporais: {len(block_time_points)}")
                print("="*80)
            if len(block_time_points) < 2: continue
            # Identificar novos elementos ativos neste bloco
            if bloco_nome == "bloco_1":
                novos_elementos = bloco_info["elementos_nos"]["elementos_dominio"]
                self._set_initial_temperature_for_new_elements(np.array(novos_elementos), self.temp_iniciais_por_camada, bloco_1=True)
                # Salvar o estado inicial (timestep 0) para conferência
                self._save_results(0, block_time_points[0])
            else:
                novos_elementos = bloco_info["diferencas"]["elementos_dominio"]["entradas"]
                if len(novos_elementos) > 0:
                    self._set_initial_temperature_for_new_elements(np.array(novos_elementos), self.temp_iniciais_por_camada, bloco_1=False)
            for i in range(1, len(block_time_points)):
                current_time = block_time_points[i]
                dt_val = current_time - block_time_points[i-1]
                total_steps += 1
                if self.rank == 0: print(f"\n[Passo {total_steps}] Tempo: {current_time/3600:.2f} h (dt = {dt_val/3600:.2f} h)")
                self._solve_timestep(dt_val, current_time)
                self._update_state()
                if i % 5 == 0 or i == len(block_time_points) - 1:
                    self._save_results(total_steps, current_time)

    def _set_initial_temperature_for_new_elements(self, novos_elementos, temp_por_pg, bloco_1=False):
        """
        Atribui a temperatura inicial correta para os elementos ativos do bloco 1 (bloco_1=True),
        ou apenas para os novos elementos ativos nos demais blocos (bloco_1=False),
        interpolando para os nós (CG-1) a partir de uma função DG-0.
        """
        temp_dg0 = Function(self.V_prop)
        # Inicializa todo o array DG-0 com zero
        temp_dg0.x.array[:] = 0.0
        if bloco_1:
            # Para o bloco 1, inicializa apenas os elementos ativos do bloco 1
            for pg_id, temp in temp_por_pg.items():
                cells = self.cell_tags.find(pg_id)
                # Só atribui para células que estão em novos_elementos (ativos no bloco 1)
                cells_ativos = np.intersect1d(cells, novos_elementos)
                temp_dg0.x.array[cells_ativos] = temp
        else:
            # Para os demais blocos, inicializa apenas os novos elementos ativos
            for pg_id, temp in temp_por_pg.items():
                cells = self.cell_tags.find(pg_id)
                cells_novos = np.intersect1d(cells, novos_elementos)
                temp_dg0.x.array[cells_novos] = temp
        # Interpola para CG-1 (nós)
        if hasattr(self, 'T'):
            self.T.interpolate(temp_dg0)
            self.T_n.interpolate(temp_dg0)
        else:
            self.Tp.interpolate(temp_dg0)
            self.Tp_n.interpolate(temp_dg0)

    def _solve_timestep(self, dt_val, current_time):
        if self.has_exothermic:
            self._update_equivalent_time_explicitly(dt_val)
            self._update_heat_generation()
        self._solve_temperature_equation(dt_val, current_time)

    def _update_equivalent_time_explicitly(self, dt_val):
        if self.rank == 0: print("   - Calculando tempo equivalente (explícito)...")
        try:
            mat_props = list(self.materials.values())[0]
            EaR, Tref = mat_props['EaR'], mat_props['Tref']
            expoente = EaR * (1 / (Tref + 273.15) - 1 / (self.Tp_n.x.array + 273.15))
            arrhenius_factor = np.exp(np.clip(expoente, -70, 70))
            self.teq.x.array[:] = self.teq_n.x.array + (dt_val * arrhenius_factor)
        except Exception as e:
            if self.rank == 0: print(f"     -> ⚠️ Erro: {e}.")

    def _update_heat_generation(self):
        if self.rank == 0: print("   - Calculando geração de calor (Q)...")
        try:
            mat_props = list(self.materials.values())[0]
            rho, ce, dTadinfty, a_sec, exp_hill = (
                mat_props['densidade'], mat_props['calor_especifico'], mat_props['dTadinfty'],
                mat_props['a_dias'] * 24 * 3600, mat_props['expoente'])
            teq_vals = self.teq.x.array
            denom = a_sec**exp_hill + teq_vals**exp_hill
            denom[denom == 0] = 1e-9
            Q_vals = rho * ce * (dTadinfty * teq_vals**exp_hill) / denom
            self.Q_heat.x.array[:] = np.maximum(Q_vals, 0.0)
        except Exception as e:
            if self.rank == 0: print(f"     -> ⚠️ Erro: {e}. Usando Q=0.")
            self.Q_heat.x.array[:] = 0.0

    def _solve_temperature_equation(self, dt_val, current_time):
        if self.rank == 0: print("   - Resolvendo equação da temperatura...")
        a, L = self._setup_variational_problem(dt_val, current_time)
        if a is None or L is None:
            if self.rank == 0: print("     -> ❌ Falha ao montar o problema variacional.")
            return
        
        # Diagnóstico adicional
        if self.rank == 0:
            print(f"     -> 🔍 Diagnóstico: dt={dt_val:.2f}s, theta={self.theta}")
            print(f"     -> 🔍 Elementos ativos: {len(self._find_active_block(current_time)['elementos_nos']['elementos_dominio'])}")
        
        bcs = self._get_boundary_conditions(current_time)
        if self.rank == 0:
            print(f"     -> 🔍 BCs aplicadas: {len(bcs)} condições de contorno")
        
        problem = LinearProblem(a, L, bcs=bcs, u=self.Tp if self.has_exothermic else self.T)
        solver_name = self._solve_with_robust_cascade(problem)
        if self.rank == 0 and solver_name:
            sol_func = self.Tp if self.has_exothermic else self.T
            T_min, T_max = np.min(sol_func.x.array), np.max(sol_func.x.array)
            print(f"     -> ✅ Convergência com '{solver_name}'. Range T: [{T_min:.2f}, {T_max:.2f}]°C")
        elif self.rank == 0:
            print(f"     -> ❌ TODOS os solvers falharam!")

    def _setup_variational_problem(self, dt_val, current_time):
        info_bloco_atual = self._find_active_block(current_time)
        if not info_bloco_atual: 
            if self.rank == 0: print("     -> ❌ Nenhum bloco ativo encontrado!")
            return None, None
        
        dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
        if self.rank == 0:
            print(f"     -> 🔍 Domínios ativos: {dominios_ativos}")
        
        dt = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        dx = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags)
        u, v = (self.u_Tp, self.v_Tp) if self.has_exothermic else (self.u, self.v)
        T_n = self.Tp_n if self.has_exothermic else self.T_n
        a, L = 0, 0
        
        for domain_id in dominios_ativos:
            a += self.rho * self.cp * u * v * dx(domain_id)
            a += dt * theta * self.k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx(domain_id)
            L += self.rho * self.cp * T_n * v * dx(domain_id)
            L -= dt * (1 - theta) * self.k * ufl.dot(ufl.grad(T_n), ufl.grad(v)) * dx(domain_id)
            if self.has_exothermic:
                L += dt * self.Q_heat * v * dx(domain_id)
        
        if self.rank == 0:
            print(f"     -> 🔍 Problema variacional montado com {len(dominios_ativos)} domínios")
        
        # Verificação adicional: garante que a e L não são zero
        if a == 0 or L == 0:
            if self.rank == 0:
                print("     -> ❌ ERRO: Forma variacional é zero!")
            return None, None
        
        return a, L

    def _get_boundary_conditions(self, current_time):
        info_bloco_atual = self._find_active_block(current_time)
        if not info_bloco_atual: return self._get_fallback_bcs()
        nos_ativos = info_bloco_atual['elementos_nos']['nos_dominio']
        contornos_ativos = info_bloco_atual['physical_groups']['lines']
        T_b = self.temp_inicial
        bcs = []
        num_total_nos = self.mesh.topology.index_map(0).size_local
        all_nodes = np.arange(num_total_nos, dtype=np.int32)
        inactive_mask = np.ones(num_total_nos, dtype=bool)
        inactive_mask[nos_ativos] = False
        inactive_dofs = all_nodes[inactive_mask]
        if inactive_dofs.size > 0:
            bcs.append(dirichletbc(Constant(self.mesh, T_b), inactive_dofs, self.V))
        fdim = self.mesh.topology.dim - 1
        for bc_id in contornos_ativos:
            facets = self.facet_tags.find(bc_id)
            if facets.size > 0:
                dofs = locate_dofs_topological(self.V, fdim, facets)
                bcs.append(dirichletbc(Constant(self.mesh, T_b), dofs, self.V))
        return bcs

    def _get_fallback_bcs(self):
        fdim = self.mesh.topology.dim - 1
        boundary_facets = mesh.exterior_facet_indices(self.mesh.topology)
        dofs = locate_dofs_topological(self.V, fdim, boundary_facets)
        return [dirichletbc(Constant(self.mesh, self.temp_inicial), dofs, self.V)]

    def _create_layer_submeshes(self):
        """
        Cria submeshes para cada camada de construção com seus respectivos physical groups.
        """
        from dolfinx.mesh import create_submesh

        self.layer_submeshes = {}
        self.layer_mappings = {}
        self.layer_boundary_mappings = {}

        # Identificar elementos por camada
        for layer_name, layer_info in self.simulation_plan.items():
            layer_id = int(layer_name.split('_')[1]) - 1  # Extrair ID da camada

            # Obter elementos ativos nesta camada
            active_elements = np.array(layer_info["elementos_nos"]["elementos_dominio"])

            if len(active_elements) > 0:
                # Criar submesh
                submesh, entity_map, vertex_map, geom_map = create_submesh(
                    self.mesh, self.mesh.topology.dim, active_elements
                )

                self.layer_submeshes[layer_id] = {
                    'submesh': submesh,
                    'entity_map': entity_map,
                    'vertex_map': vertex_map,
                    'geom_map': geom_map,
                    'active_elements': active_elements
                }

                # Mapear physical groups para esta submesh
                self._map_boundary_groups_to_submesh(layer_id, layer_info)

    def _map_boundary_groups_to_submesh(self, layer_id, layer_info):
        """
        Mapeia os physical groups (contornos) para cada submesh.
        """
        active_contours = layer_info["physical_groups"]["lines"]
        submesh = self.layer_submeshes[layer_id]['submesh']

        # Criar facet tags para a submesh
        fdim = submesh.topology.dim - 1
        submesh.topology.create_entities(fdim)
        submesh.topology.create_connectivity(fdim, submesh.topology.dim)

        # Mapear facets da submesh para facets da malha original
        boundary_mappings = {}

        for contour_id in active_contours:
            # Encontrar facets na malha original
            original_facets = self.facet_tags.find(contour_id)

            # Mapear para facets na submesh
            if len(original_facets) > 0:
                # Identificar quais facets pertencem à submesh
                submesh_facets = []
                for facet_idx in original_facets:
                    # Verificar se a facet está na fronteira da submesh
                    # Isso requer lógica mais sofisticada para mapeamento correto
                    submesh_facets.append(facet_idx)

                boundary_mappings[contour_id] = {
                    'original_facets': original_facets,
                    'submesh_facets': np.array(submesh_facets, dtype=np.int32)
                }

        self.layer_boundary_mappings[layer_id] = boundary_mappings

    def _get_boundary_conditions_for_submesh(self, layer_id, current_time):
        """
        Aplica condições de contorno robustas para uma submesh específica.

        Args:
            layer_id: ID da camada/submesh
            current_time: Tempo atual da simulação

        Returns:
            Lista de condições de contorno para a submesh
        """
        if layer_id not in self.layer_submeshes:
            return []

        submesh_info = self.layer_submeshes[layer_id]
        submesh = submesh_info['submesh']

        # Criar espaço de funções para a submesh
        V_sub = functionspace(submesh, ("Lagrange", 1))

        # Obter configurações de contorno para esta camada
        layer_name = f"camada_{layer_id + 1}"
        active_contours = self._get_active_contours_for_layer(layer_name, current_time)

        bcs = []

        # Processar cada contorno ativo
        for contour_config in active_contours:
            contour_id = contour_config['id']
            contour_type = contour_config['tipo']

            # Obter facets na submesh
            if contour_id in self.layer_boundary_mappings[layer_id]:
                submesh_facets = self.layer_boundary_mappings[layer_id][contour_id]['submesh_facets']

                if len(submesh_facets) > 0:
                    # Calcular valor da condição de contorno
                    bc_value = self._calculate_boundary_value(contour_config, current_time)

                    # Aplicar ao elemento (DG-0) e interpolar para nós (CG-1)
                    bc_function = self._apply_boundary_to_elements_then_interpolate(
                        submesh, submesh_facets, bc_value, contour_type, V_sub
                    )

                    # Criar condição de contorno nos nós
                    fdim = submesh.topology.dim - 1
                    dofs = locate_dofs_topological(V_sub, fdim, submesh_facets)
                    bcs.append(dirichletbc(bc_function, dofs, V_sub))

        return bcs

    def _get_active_contours_for_layer(self, layer_name, current_time):
        """
        Identifica quais contornos estão ativos para uma camada específica.
        """
        active_contours = []

        for contour in self.config.get('contornos', []):
            # Verificar se o contorno nasce com esta camada
            if contour.get('nasce_com_camada') == layer_name:
                # Verificar se está desativado
                desativado_em = contour.get('desativado_em', None)
                if desativado_em is None or current_time < desativado_em:
                    active_contours.append(contour)

        return active_contours

    def _calculate_boundary_value(self, contour_config, current_time):
        """
        Calcula o valor da condição de contorno baseado no tipo e parâmetros.
        """
        contour_type = contour_config['tipo']

        if contour_type == "dirichlet":
            return contour_config.get('valor', 20.0)

        elif contour_type == "conveccao":
            h = contour_config.get('h', 8.0)
            t_ext = contour_config.get('t_ext', 25.0)
            # Para convecção, retornamos t_ext como condição de Dirichlet equivalente
            return t_ext

        elif contour_type == "fluxo":
            flux_value = contour_config.get('fluxo', 0.0)
            # Para fluxo, retornamos temperatura ambiente como aproximação
            return contour_config.get('t_ambiente', 20.0)

        elif contour_type == "robin":
            h = contour_config.get('h', 8.0)
            t_ext = contour_config.get('t_ext', 25.0)
            return t_ext

        else:
            return 20.0  # Valor padrão

    def _apply_boundary_to_elements_then_interpolate(self, submesh, facets, bc_value,
                                                   contour_type, V_sub):
        """
        Aplica condição de contorno aos elementos e interpola para os nós.

        Args:
            submesh: Submesh onde aplicar a condição
            facets: Facets do contorno
            bc_value: Valor da condição de contorno
            contour_type: Tipo de condição (dirichlet, conveccao, etc.)
            V_sub: Espaço de funções da submesh

        Returns:
            Função com valores interpolados nos nós
        """
        # Criar função DG-0 para valores por elemento
        V_DG = functionspace(submesh, ("DG", 0))
        bc_dg = Function(V_DG)

        # Identificar elementos adjacentes às facets do contorno
        boundary_elements = self._get_boundary_elements(submesh, facets)

        # Atribuir valor aos elementos do contorno
        if len(boundary_elements) > 0:
            bc_dg.x.array[boundary_elements] = bc_value

        # Interpolar para espaço CG-1 (nós)
        bc_cg = Function(V_sub)
        bc_cg.interpolate(bc_dg)

        return bc_cg

    def _get_boundary_elements(self, submesh, boundary_facets):
        """
        Identifica elementos que têm facets no contorno.
        """
        # Obter conectividade entre facets e células
        tdim = submesh.topology.dim
        fdim = tdim - 1

        # Criar conectividade se necessário
        submesh.topology.create_connectivity(fdim, tdim)
        facet_to_cell = submesh.topology.connectivity(fdim, tdim)

        boundary_elements = []

        for facet_idx in boundary_facets:
            if facet_idx < len(facet_to_cell):
                cells = facet_to_cell.links(facet_idx)
                boundary_elements.extend(cells)

        return np.unique(np.array(boundary_elements, dtype=np.int32))

    def _synchronize_boundary_conditions(self, layer_id, bcs_submesh):
        """
        Sincroniza condições de contorno entre submesh e malha principal.
        """
        if layer_id not in self.layer_mappings:
            return

        # Mapear de volta para a malha principal se necessário
        # Implementação específica depende da necessidade de sincronização
        pass

    def _solve_with_robust_cascade(self, problem):
        # Configurações mais conservadoras para evitar problemas de memória
        solvers = {
            "GMRES+ILU (conservador)": {
                "ksp_type": "gmres", 
                "pc_type": "ilu", 
                "ksp_rtol": 1e-6,
                "ksp_max_it": 100,
                "ksp_gmres_restart": 30
            },
            "CG+JACOBI": {
                "ksp_type": "cg", 
                "pc_type": "jacobi", 
                "ksp_rtol": 1e-6,
                "ksp_max_it": 200
            },
            "LU (direto)": {
                "ksp_type": "preonly", 
                "pc_type": "lu"
            }
        }
        
        for name, opts in solvers.items():
            try:
                # Limpa opções anteriores
                problem.petsc_options = {}
                # Aplica novas opções
                problem.petsc_options = opts
                problem.solve()
                return name
            except (RuntimeError, PETSc.Error) as e:
                if self.rank == 0: 
                    print(f"     -> ⚠️ Solver {name} falhou: {str(e)[:100]}...")
            except Exception as e:
                if self.rank == 0:
                    print(f"     -> ⚠️ Solver {name} falhou com erro inesperado: {str(e)[:100]}...")
        
        # Se todos falharam, tenta uma abordagem mais simples
        try:
            if self.rank == 0: print("     -> 🔄 Tentando abordagem simplificada...")
            problem.petsc_options = {
                "ksp_type": "gmres",
                "pc_type": "none",
                "ksp_rtol": 1e-4,
                "ksp_max_it": 50
            }
            problem.solve()
            return "GMRES (simplificado)"
        except Exception as e:
            if self.rank == 0:
                print(f"     -> ❌ Abordagem simplificada também falhou: {str(e)[:100]}...")
        
        return None

    def _update_state(self):
        if self.has_exothermic:
            self.Tp_n.x.array[:] = self.Tp.x.array
            self.teq_n.x.array[:] = self.teq.x.array
        else:
            self.T_n.x.array[:] = self.T.x.array

    def _find_active_block(self, current_time):
        for block_data in self.simulation_plan.values():
            if block_data['info_bloco']['inicio'] <= current_time <= block_data['info_bloco']['fim']:
                return block_data
        return None

    def _save_results(self, time_step_idx, current_time):
        """Salva resultados principais e também submeshes se ativadas."""
        if self.rank == 0:
            print(f"   💾 Salvando passo {time_step_idx} (t={current_time/3600:.2f}h)...")
        try:
            funcs = [self.Tp, self.teq, self.Q_heat] if self.has_exothermic else [self.T]
            for func in funcs:
                filename = Path(self.output_dir) / f"{func.name}_passo_{time_step_idx:04d}.xdmf"
                with io.XDMFFile(self.comm, filename, "w") as xdmf:
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(func, current_time)

            # Salvar submeshes e resultados por camada
            if hasattr(self, 'layer_submeshes') and self.layer_submeshes:
                self._save_all_submeshes_and_results(time_step_idx, current_time)

        except Exception as e:
            if self.rank == 0: print(f"     -> ⚠️ Erro ao salvar: {e}")
            
    def _finalize(self):
        """Finaliza a simulação e salva informações adicionais."""
        if self.rank == 0:
            print("\n" + "="*80)
            print("✅ SIMULAÇÃO CONCLUÍDA")
            print(f"   - Resultados salvos em: '{self.output_dir}'")

        # Salvar submeshes na pasta do caso
        case_dir = Path(self.config_file).parent
        self._save_submeshes(case_dir)

        if self.rank == 0:
            print("   ✅ Submeshes salvas na pasta do caso.")
            print("   ✅ Simulação concluída com sucesso!")
            print(f"   📁 Resultados salvos em: {case_dir}")

    def _save_submeshes(self, case_dir):
        """
        Salva todas as submeshes em arquivos XDMF separados na pasta do caso.

        Args:
            case_dir: Diretório do caso onde salvar as submeshes
        """
        if not hasattr(self, 'layer_submeshes'):
            return

        if self.rank == 0:
            print(f"\n--- SALVANDO SUBMESHES ---")
            print(f"   📁 Diretório: {case_dir}")

        submesh_dir = case_dir / "submeshes"
        submesh_dir.mkdir(exist_ok=True)

        for layer_id, submesh_info in self.layer_submeshes.items():
            submesh = submesh_info['submesh']
            layer_name = f"camada_{layer_id + 1}"

            # Nome do arquivo baseado no caso
            caso_name = case_dir.name
            filename = submesh_dir / f"{caso_name}_{layer_name}.xdmf"

            if self.rank == 0:
                print(f"   ➡️  Salvando {layer_name}: {filename.name}")

            # Salvar submesh
            with io.XDMFFile(self.comm, str(filename), "w") as xdmf:
                xdmf.write_mesh(submesh)

                # Salvar cell tags se existirem
                if hasattr(self, 'cell_tags') and self.cell_tags is not None:
                    # Filtrar cell tags para esta submesh
                    active_elements = submesh_info['active_elements']
                    filtered_tags = self._filter_cell_tags_for_submesh(
                        self.cell_tags, active_elements, submesh_info['entity_map']
                    )
                    if filtered_tags is not None:
                        xdmf.write_meshtags(filtered_tags, submesh.geometry)

                # Salvar facet tags (contornos) para esta submesh
                boundary_tags = self._create_boundary_tags_for_submesh(layer_id)
                if boundary_tags is not None:
                    xdmf.write_meshtags(boundary_tags, submesh.geometry)

            if self.rank == 0:
                print(f"   ✅ {layer_name} salva com sucesso!")

    def _filter_cell_tags_for_submesh(self, original_tags, active_elements, entity_map):
        """
        Filtra cell tags para incluir apenas elementos da submesh.

        Args:
            original_tags: Cell tags da malha original
            active_elements: Elementos ativos na submesh
            entity_map: Mapeamento de entidades da submesh para a malha original

        Returns:
            Cell tags filtradas para a submesh ou None se não aplicável
        """
        if original_tags is None:
            return None

        try:
            # Criar novo array de tags para a submesh
            submesh_tags = np.full(len(entity_map), -1, dtype=np.int32)

            # Mapear tags da malha original para a submesh
            for sub_idx, orig_idx in enumerate(entity_map):
                if orig_idx < len(original_tags.values):
                    submesh_tags[sub_idx] = original_tags.values[orig_idx]

            # Criar meshtags para a submesh
            from dolfinx.mesh import meshtags
            submesh = self.layer_submeshes[list(self.layer_submeshes.keys())[0]]['submesh']

            # Criar tags válidas (remover -1)
            valid_mask = submesh_tags != -1
            if not np.any(valid_mask):
                return None

            # Criar meshtags com as tags válidas
            return meshtags(submesh, submesh.topology.dim,
                          np.arange(len(submesh_tags), dtype=np.int32),
                          submesh_tags)

        except Exception as e:
            if self.rank == 0:
                print(f"   ⚠️  Erro ao filtrar cell tags: {e}")
            return None

    def _create_boundary_tags_for_submesh(self, layer_id):
        """
        Cria facet tags (contornos) para uma submesh específica.

        Args:
            layer_id: ID da camada

        Returns:
            Facet tags para a submesh ou None se não aplicável
        """
        if layer_id not in self.layer_boundary_mappings:
            return None

        try:
            submesh = self.layer_submeshes[layer_id]['submesh']
            boundary_mappings = self.layer_boundary_mappings[layer_id]

            # Coletar todas as facets e suas tags
            all_facets = []
            all_tags = []

            for contour_id, mapping in boundary_mappings.items():
                submesh_facets = mapping.get('submesh_facets', [])
                if len(submesh_facets) > 0:
                    all_facets.extend(submesh_facets)
                    all_tags.extend([contour_id] * len(submesh_facets))

            if not all_facets:
                return None

            # Criar arrays numpy
            facets_array = np.array(all_facets, dtype=np.int32)
            tags_array = np.array(all_tags, dtype=np.int32)

            # Criar meshtags
            from dolfinx.mesh import meshtags
            fdim = submesh.topology.dim - 1
            return meshtags(submesh, fdim, facets_array, tags_array)

        except Exception as e:
            if self.rank == 0:
                print(f"   ⚠️  Erro ao criar boundary tags: {e}")
            return None

    def _save_submesh_results(self, layer_id, solution, step, time):
        """
        Salva resultados de uma submesh específica.

        Args:
            layer_id: ID da camada
            solution: Função de solução
            step: Número do passo
            time: Tempo atual
        """
        case_dir = Path(self.config_file).parent
        submesh_dir = case_dir / "submeshes" / "results"
        submesh_dir.mkdir(parents=True, exist_ok=True)

        layer_name = f"camada_{layer_id + 1}"
        caso_name = case_dir.name

        filename = submesh_dir / f"{caso_name}_{layer_name}_step{step:04d}.xdmf"

        submesh = self.layer_submeshes[layer_id]['submesh']

        with io.XDMFFile(self.comm, str(filename), "w") as xdmf:
            xdmf.write_mesh(submesh)
            xdmf.write_function(solution, time)

    def _save_all_submeshes_and_results(self, step, current_time):
        """
        Função auxiliar para salvar todas as submeshes e resultados atuais.

        Args:
            step: Número do passo
            current_time: Tempo atual da simulação
        """
        case_dir = Path(self.config_file).parent

        # Salvar estrutura das submeshes (apenas uma vez)
        if step == 0:
            self._save_submeshes(case_dir)

        # Salvar resultados de cada submesh ativa
        for layer_id in self.layer_submeshes.keys():
            if hasattr(self, 'T_layers') and layer_id in self.T_layers:
                self._save_submesh_results(layer_id, self.T_layers[layer_id], step, current_time)
            elif hasattr(self, 'T') and self.T is not None:
                # Se não houver solução separada por camada, usar a principal
                self._save_submesh_results(layer_id, self.T, step, current_time)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"USO: python {sys.argv[0]} <caminho_para_pasta_do_caso>")
        sys.exit(1)
    pasta_caso = Path(sys.argv[1])
    caso = pasta_caso.name
    config_file = pasta_caso / f"{caso}.yaml"
    json_file = pasta_caso / f"{caso}-xdmf.json"
    log_file = pasta_caso / "log_simulacao.md"
    sim = SimulacaoBarragemR2(str(config_file), str(json_file), str(log_file))
    sim.run()