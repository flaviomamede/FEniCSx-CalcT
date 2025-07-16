#!/usr/bin/env python3
"""
BARRAGEM FEniCSx - VERSÃO FINAL REATORADA

Este script implementa uma análise de transferência de calor em uma barragem
construída em etapas (stagewise), utilizando uma abordagem robusta e modular.

PRINCIPAIS CARACTERÍSTICAS:
- **Estrutura Orientada a Objetos:** A lógica da simulação é encapsulada na
  classe `SimulacaoBarragem` para melhor organização.
- **Acionamento por JSON:** A simulação é inteiramente guiada pelo arquivo
  `analise_stagewise_xdmf.json`, que define os domínios, contornos e
  camadas ativas em cada bloco de tempo.
- **Lógica de Solver Robusta:** Implementa uma cascata de solvers (GMRES -> CG -> LU)
  para garantir a convergência mesmo em cenários desafiadores.
- **Tratamento de Domínio em Crescimento:** Utiliza a técnica de fixar os graus
  de liberdade (DOFs) de nós inativos para evitar matrizes singulares, uma
  abordagem essencial para problemas de construção em etapas.
- **Clareza e Manutenibilidade:** Métodos com responsabilidades únicas e
  docstrings detalhadas para facilitar o entendimento, a depuração e futuros
  desenvolvimentos.
- **Extensibilidade:** Preparado para executar a simulação completa em todos os
  blocos de tempo com uma modificação mínima.
"""

import os
import sys
import json
import yaml
import datetime
import numpy as np
import ufl
import dolfinx as fem
from dolfinx import mesh, io
from dolfinx.fem import Function, FunctionSpace, Constant, locate_dofs_topological
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

class SimulacaoBarragem:
    """
    Encapsula toda a lógica para a simulação de transferência de calor
    em uma barragem construída em etapas.
    """

    def __init__(self, config_file, json_file, log_file="log_simulacao.md"):
        """
        Inicializa a simulação, carregando configurações e preparando o ambiente.

        Args:
            config_file (str): Caminho para o arquivo de configuração .yaml.
            json_file (str): Caminho para o arquivo de análise .json.
            log_file (str): Nome do arquivo de log a ser gerado.
        """
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        
        # Redireciona a saída para o log
        self.logger = TeeOutput(log_file)
        sys.stdout = self.logger

        if self.rank == 0:
            print("="*80)
            print("🏗️  INICIALIZANDO SIMULAÇÃO DE BARRAGEM STAGE-WISE COM FENICSx")
            print("="*80)

        self._load_config(config_file, json_file)

    def _load_config(self, config_file, json_file):
        """Carrega os arquivos de configuração YAML e JSON."""
        if self.rank == 0:
            print(f"🔄 Carregando configuração de '{config_file}' e análise de '{json_file}'...")

        with open(config_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        with open(json_file, 'r', encoding='utf-8') as f:
            self.analysis = json.load(f)

        # Parâmetros gerais
        self.theta = self.config['general']['theta']
        self.output_dir = self.config['general']['output_dir']
        
        # Ajustar o caminho da malha para o diretório do caso
        mesh_file_rel = self.config['general']['mesh_file']
        case_dir = Path(config_file).parent
        self.mesh_file_path = case_dir / mesh_file_rel

        # Informações da análise JSON
        self.time_vector = self.analysis['vetor_tempo']
        self.time_blocks = self.analysis['blocos_tempo']
        self.simulation_plan = self.analysis['analise_resultados']
        self.mappings = self.analysis['info_geral']['mapeamentos']

        # Temperatura inicial (lógica de fallback)
        self.temp_inicial = 20.0
        if 'initial_conditions' in self.config and 'temperature' in self.config['initial_conditions']:
            self.temp_inicial = self.config['initial_conditions']['temperature']
        elif 'camadas_material' in self.config and self.config['camadas_material']:
            # Pega a temperatura da primeira camada como referência
            self.temp_inicial = self.config['camadas_material'][0]['temperatura_inicial']

        if self.rank == 0:
            print(f"   ✅ Configuração carregada com sucesso.")
            print(f"   🌡️ Temperatura Inicial Padrão: {self.temp_inicial}°C")
            print(f"   📅 Número de Blocos Construtivos: {len(self.time_blocks)}")
            print(f"   📁 Diretório de Saída: {self.output_dir}")
            print(f"   📁 Arquivo de Malha: {self.mesh_file_path}")

    def run(self):
        """
        Orquestra a execução completa da simulação.
        """
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
            if self.logger:
                self.logger.close()

    def _setup(self):
        """Executa todas as rotinas de configuração inicial."""
        if self.rank == 0:
            print("\n--- FASE DE CONFIGURAÇÃO ---")
        self._load_mesh()
        self._setup_materials()
        self._setup_function_spaces()
        self._set_initial_conditions()
        os.makedirs(self.output_dir, exist_ok=True)

    def _load_mesh(self):
        """Carrega a malha e as tags de subdomínio (células) e contorno (facetas)."""
        if self.rank == 0:
            print(f"   ➡️  Carregando malha de '{self.mesh_file_path}'...")
        
        with io.XDMFFile(self.comm, self.mesh_file_path, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            self.mesh.topology.create_entities(1) # Garante a criação de facetas
            self.mesh.topology.create_connectivity(1, 2) # Conectividade faceta -> célula
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        
        if self.rank == 0:
            print(f"   ✅ Malha '{self.mesh.name}' carregada.")
            print(f"      - Células: {self.mesh.topology.index_map(self.mesh.topology.dim).size_global}")
            print(f"      - Vértices: {self.mesh.topology.index_map(0).size_global}")
            print(f"      - Tags de Célula (Volumes/Surfaces) detectadas.")
            print(f"      - Tags de Faceta (Surfaces/Lines) detectadas.")

    def _setup_materials(self):
        """Configura as propriedades dos materiais a partir do arquivo YAML."""
        if self.rank == 0:
            print("   ➡️  Configurando propriedades dos materiais...")
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
        if self.rank == 0:
            print(f"   ✅ Materiais configurados: {list(self.materials.keys())}")
            if self.has_exothermic:
                print("   🔥 Reação exotérmica detectada. O sistema acoplado (Temperatura + Tempo Equivalente) será utilizado.")
            else:
                print("   ❄️ Apenas transferência de calor simples será considerada.")

    def _setup_function_spaces(self):
        """Define os espaços de função necessários para o problema."""
        if self.rank == 0:
            print("   ➡️  Definindo espaços de função (Lagrange P1)...")
        
        # Corrigir a sintaxe do FunctionSpace para a versão atual do FEniCSx
        self.V = FunctionSpace(self.mesh, ("CG", 1))

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

    def _set_initial_conditions(self):
        """Aplica as condições iniciais em todas as funções."""
        if self.rank == 0:
            print("   ➡️  Aplicando condições iniciais...")
        
        if self.has_exothermic:
            self.Tp_n.x.array[:] = self.temp_inicial
            self.Tp.x.array[:] = self.temp_inicial
            self.teq_n.x.array[:] = 0.0
            self.teq.x.array[:] = 0.0
            self.Q_heat.x.array[:] = 0.0
            if self.rank == 0:
                print(f"   ✅ Condições iniciais (exotérmico): T={self.temp_inicial}°C, teq=0s")
        else:
            self.T_n.x.array[:] = self.temp_inicial
            self.T.x.array[:] = self.temp_inicial
            if self.rank == 0:
                print(f"   ✅ Condições iniciais (simples): T={self.temp_inicial}°C")

    def _run_simulation_loop(self):
        """
        Executa o loop principal da simulação através dos blocos de tempo e passos de tempo.
        """
        if self.rank == 0:
            print("\n--- FASE DE SIMULAÇÃO ---")
            print("🚀 Iniciando loop temporal...")
        
        total_steps = 0
        # NOTA: Para rodar apenas o primeiro bloco, use: self.time_blocks[:1]
        # Para rodar a simulação completa, use: self.time_blocks
        for block_idx, block_info in enumerate(self.time_blocks[:1]):
            
            # Filtra o vetor de tempo global para os pontos relevantes deste bloco
            block_time_points = [t for t in self.time_vector if block_info['inicio'] <= t <= block_info['fim']]
            
            if self.rank == 0:
                print("\n" + "="*80)
                print(f"📦 PROCESSANDO BLOCO CONSTRUTIVO {block_idx + 1}")
                print(f"   - Período: {block_info['inicio']/3600:.1f}h a {block_info['fim']/3600:.1f}h")
                print(f"   - Duração: {block_info['duracao']/3600:.1f}h")
                print(f"   - Passos de Tempo no Bloco: {len(block_time_points) - 1}")
                print("="*80)

            if len(block_time_points) < 2:
                continue

            for i in range(1, len(block_time_points)):
                current_time = block_time_points[i]
                previous_time = block_time_points[i-1]
                dt_val = current_time - previous_time
                total_steps += 1

                if self.rank == 0:
                    print(f"\n[Passo {total_steps}] Tempo: {current_time/3600:.2f} h (dt = {dt_val/3600:.2f} h)")
                
                self._solve_timestep(dt_val, current_time)
                self._update_state()
                
                # Salvar resultados periodicamente
                if i % 5 == 0 or i == len(block_time_points) - 1:
                    self._save_results(total_steps, current_time)

    def _solve_timestep(self, dt_val, current_time):
        """
        Resolve um único passo de tempo.

        Args:
            dt_val (float): O incremento de tempo (delta t).
            current_time (float): O tempo absoluto atual da simulação.
        """
        # Abordagem de Staggered/Operator-Splitting para o problema exotérmico
        if self.has_exothermic:
            self._update_equivalent_time_explicitly(dt_val, current_time)
            self._update_heat_generation(current_time)
        
        # Resolve a equação da temperatura (com ou sem geração de calor)
        self._solve_temperature_equation(dt_val, current_time)

    def _update_equivalent_time_explicitly(self, dt_val, current_time):
        """
        PASSO 1 (Exotérmico): Atualiza o tempo equivalente de forma explícita,
        baseado na temperatura do passo de tempo anterior (Tp_n).
        teq(t+dt) = teq(t) + dt * Arrhenius(Tp(t))
        """
        if self.rank == 0:
            print("   - Calculando tempo equivalente (explícito)...")
        
        try:
            # Propriedades do primeiro material (simplificação, idealmente seria por domínio)
            mat_props = list(self.materials.values())[0]
            EaR, Tref = mat_props['EaR'], mat_props['Tref']
            
            Tp_anterior_vals = self.Tp_n.x.array
            
            # Fator de Arrhenius: exp(EaR * (1/(Tref_K) - 1/(Tp_K)))
            expoente = EaR * (1 / (Tref + 273.15) - 1 / (Tp_anterior_vals + 273.15))
            expoente_clipped = np.clip(expoente, -70, 70) # Previne overflow/underflow em exp()
            arrhenius_factor = np.exp(expoente_clipped)
            
            # Atualização explícita
            teq_increment = dt_val * arrhenius_factor
            self.teq.x.array[:] = self.teq_n.x.array + teq_increment

            if self.rank == 0:
                teq_min, teq_max = np.min(self.teq.x.array), np.max(self.teq.x.array)
                print(f"     -> Tempo equivalente atualizado. Range: [{teq_min/3600:.1f}h, {teq_max/3600:.1f}h]")

        except Exception as e:
            if self.rank == 0:
                print(f"     -> ⚠️ Erro ao atualizar tempo equivalente: {e}. Mantendo valor anterior.")
            self.teq.x.array[:] = self.teq_n.x.array

    def _update_heat_generation(self, current_time):
        """
        PASSO 2 (Exotérmico): Calcula o campo de geração de calor (Q_heat)
        baseado no tempo equivalente recém-calculado.
        """
        if self.rank == 0:
            print("   - Calculando geração de calor (Q)...")
            
        try:
            mat_props = list(self.materials.values())[0]
            rho, ce = mat_props['densidade'], mat_props['calor_especifico']
            dTadinfty = mat_props['dTadinfty']
            a_sec = mat_props['a_dias'] * 24 * 3600
            expoente = mat_props['expoente']

            teq_vals = self.teq.x.array
            
            # Função de geração de calor (Hill)
            numerador = dTadinfty * (teq_vals**expoente)
            denominador = a_sec**expoente + teq_vals**expoente
            
            # Evita divisão por zero se teq_vals for negativo (não físico)
            denominador[denominador == 0] = 1e-9

            Q_vals = rho * ce * (numerador / denominador)
            self.Q_heat.x.array[:] = np.maximum(Q_vals, 0.0) # Garante não-negatividade

            if self.rank == 0:
                q_min, q_max = np.min(self.Q_heat.x.array), np.max(self.Q_heat.x.array)
                print(f"     -> Geração de calor calculada. Range: [{q_min:.1e}, {q_max:.1e}] W/m³")

        except Exception as e:
            if self.rank == 0:
                print(f"     -> ⚠️ Erro ao calcular geração de calor: {e}. Usando Q=0.")
            self.Q_heat.x.array[:] = 0.0

    def _solve_temperature_equation(self, dt_val, current_time):
        """
        PASSO 3: Monta e resolve o sistema linear para a equação da temperatura.
        """
        if self.rank == 0:
            print("   - Resolvendo equação da temperatura...")
            
        # 1. Montar o problema variacional (lados esquerdo 'a' e direito 'L')
        a, L = self._setup_variational_problem(dt_val, current_time)
        if a is None:
            if self.rank == 0:
                print("     -> ❌ Falha ao montar o problema variacional. Abortando passo.")
            # Se falhar, a temperatura atual (T ou Tp) não será alterada,
            # e no próximo passo, T_n será igual ao T_n anterior.
            return

        # 2. Obter condições de contorno
        bcs = self._get_boundary_conditions(current_time)

        # 3. Resolver o sistema linear
        problem = LinearProblem(a, L, bcs=bcs, u=self.Tp if self.has_exothermic else self.T)
        solver_name = self._solve_with_robust_cascade(problem)
        
        # 4. Relatar resultado
        if self.rank == 0:
            if solver_name:
                sol_func = self.Tp if self.has_exothermic else self.T
                T_min, T_max = np.min(sol_func.x.array), np.max(sol_func.x.array)
                print(f"     -> ✅ Convergência alcançada com solver '{solver_name}'.")
                print(f"     -> 🌡️ Range de Temperatura: [{T_min:.2f}, {T_max:.2f}]°C")
            else:
                print("     -> ❌ FALHA: Nenhum solver conseguiu convergir. A solução do passo anterior será mantida.")

    def _setup_variational_problem(self, dt_val, current_time):
        """
        Monta a forma variacional (bilinear 'a' e linear 'L') para a equação de calor.
        Usa a lógica do JSON para ativar apenas os domínios relevantes.
        """
        info_bloco_atual = self._find_active_block(current_time)
        if not info_bloco_atual:
            if self.rank == 0:
                print(f"     -> ⚠️ Nenhum bloco de simulação ativo para t={current_time/3600:.1f}h. Pulando montagem.")
            return None, None

        dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
        if self.rank == 0:
            print(f"     -> Domínios ativos para montagem: {dominios_ativos}")

        dt = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        
        # Simplificação: usa o primeiro material. Idealmente, as propriedades
        # (k, rho, cp) deveriam ser Funções Discretas (DG-0) definidas por domínio.
        mat = list(self.materials.values())[0]
        k = Constant(self.mesh, PETSc.ScalarType(mat['condutividade']))
        rho = Constant(self.mesh, PETSc.ScalarType(mat['densidade']))
        cp = Constant(self.mesh, PETSc.ScalarType(mat['calor_especifico']))
        
        # Medida de integração sobre os subdomínios (células)
        dx_measure = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags)

        # Seleciona as funções corretas (exotérmico vs. simples)
        u, v = (self.u_Tp, self.v_Tp) if self.has_exothermic else (self.u, self.v)
        T_n = self.Tp_n if self.has_exothermic else self.T_n
        
        # Inicializa formas
        a, L = 0, 0
        
        # Contribuição de cada domínio ativo
        for domain_id in dominios_ativos:
            # Lado Esquerdo (Bilinear form 'a')
            a += rho * cp * u * v * dx_measure(domain_id)
            a += dt * theta * k * ufl.dot(ufl.grad(u), ufl.grad(v)) * dx_measure(domain_id)

            # Lado Direito (Linear form 'L')
            L += rho * cp * T_n * v * dx_measure(domain_id)
            L -= dt * (1 - theta) * k * ufl.dot(ufl.grad(T_n), ufl.grad(v)) * dx_measure(domain_id)
            
            # Termo fonte de calor (apenas para sistema exotérmico)
            if self.has_exothermic:
                L += dt * self.Q_heat * v * dx_measure(domain_id)
        
        # Retorna o lado esquerdo e direito completos da equação: a = dt*L
        return a, L

    def _get_boundary_conditions(self, current_time):
        """
        Define as condições de contorno de Dirichlet. Esta é a parte mais crítica
        para a análise stagewise.

        ESTRATÉGIA:
        1. Identifica todos os nós *inativos* (fora do domínio de simulação atual)
           e aplica uma condição de Dirichlet neles para fixar seu valor,
           evitando uma matriz singular.
        2. Aplica as condições de Dirichlet nos contornos *ativos* definidos no JSON.
        """
        if self.rank == 0:
            print("   - Aplicando Condições de Contorno...")
            
        info_bloco_atual = self._find_active_block(current_time)
        if not info_bloco_atual:
            if self.rank == 0: print("     -> ⚠️ Usando BCs de fallback (fronteira completa).")
            return self._get_fallback_bcs()

        nos_ativos = info_bloco_atual['elementos_nos']['nos_dominio']
        contornos_ativos_ids = info_bloco_atual['physical_groups']['lines']
        
        T_boundary = self.temp_inicial # Simplificação: usa a temperatura inicial
        
        bcs = []
        
        # 1. Fixar nós inativos (essencial para estabilidade)
        num_total_nos = self.mesh.topology.index_map(0).size_local
        all_nodes = np.arange(num_total_nos, dtype=np.int32)
        inactive_nodes_mask = np.ones(num_total_nos, dtype=bool)
        inactive_nodes_mask[nos_ativos] = False
        inactive_nodes = all_nodes[inactive_nodes_mask]

        if inactive_nodes.size > 0:
            # Os nós são os DOFs para o espaço Lagrange P1
            inactive_dofs = inactive_nodes
            bc_inactive = fem.dirichletbc(Constant(self.mesh, T_boundary), inactive_dofs, self.V)
            bcs.append(bc_inactive)
            if self.rank == 0: print(f"     -> {len(inactive_dofs)} DOFs inativos fixados em T={T_boundary}°C.")

        # 2. Aplicar BCs nos contornos ativos
        if self.rank == 0: print(f"     -> Contornos ativos para BCs: {contornos_ativos_ids}")
        for boundary_id in contornos_ativos_ids:
            facet_indices = np.where(self.facet_tags.values == boundary_id)[0]
            if facet_indices.size > 0:
                boundary_dofs = locate_dofs_topological(self.V, self.mesh.topology.dim - 1, facet_indices)
                bc = fem.dirichletbc(Constant(self.mesh, T_boundary), boundary_dofs, self.V)
                bcs.append(bc)
        
        if self.rank == 0: print(f"     -> Total de {len(bcs)} objetos de BC criados.")
        return bcs

    def _get_fallback_bcs(self):
        """Retorna BCs para toda a fronteira externa, caso a lógica principal falhe."""
        boundary_facets = mesh.exterior_facet_indices(self.mesh.topology)
        boundary_dofs = locate_dofs_topological(self.V, self.mesh.topology.dim - 1, boundary_facets)
        return [fem.dirichletbc(Constant(self.mesh, self.temp_inicial), boundary_dofs, self.V)]

    def _solve_with_robust_cascade(self, problem):
        """
        Tenta resolver o sistema linear usando uma sequência de solvers,
        do mais rápido/menos robusto para o mais lento/mais robusto.

        Returns:
            str: O nome do solver que obteve sucesso, ou None se todos falharem.
        """
        solvers = {
            "GMRES+ILU": {"ksp_type": "gmres", "pc_type": "ilu", "ksp_rtol": 1e-8, "ksp_atol": 1e-10},
            "CG+HYPRE": {"ksp_type": "cg", "pc_type": "hypre", "ksp_rtol": 1e-8, "ksp_atol": 1e-10},
            "LU (direto)": {"ksp_type": "preonly", "pc_type": "lu", "pc_factor_mat_solver_type": "mumps"}
        }

        for name, opts in solvers.items():
            try:
                if self.rank == 0: print(f"     -> Tentando solver: {name}...")
                problem.petsc_options = opts
                solution = problem.solve()
                # Checa por NaN ou Inf na solução
                if not np.all(np.isfinite(solution.x.array)):
                    raise RuntimeError("Solução contém valores não finitos (NaN/Inf).")
                return name
            except (RuntimeError, PETSc.Error) as e:
                if self.rank == 0: print(f"     -> ⚠️ Solver {name} falhou: {str(e).splitlines()[0]}")
        
        return None

    def _update_state(self):
        """
        Atualiza as funções do passo anterior (_n) com a solução recém-calculada
        para preparar para o próximo passo de tempo.
        """
        if self.has_exothermic:
            self.Tp_n.x.array[:] = self.Tp.x.array
            self.teq_n.x.array[:] = self.teq.x.array
        else:
            self.T_n.x.array[:] = self.T.x.array

    def _find_active_block(self, current_time):
        """
        Encontra o bloco de construção ativo no `simulation_plan` para um dado tempo.
        """
        for block_data in self.simulation_plan.values():
            info = block_data['info_bloco']
            if info['inicio'] <= current_time <= info['fim']:
                return block_data
        return None

    def _save_results(self, time_step_idx, current_time):
        """Salva os campos de resultado em formato XDMF."""
        if self.rank == 0:
            print(f"   💾 Salvando resultados para o passo {time_step_idx} (t={current_time/3600:.2f}h)...")
        try:
            if self.has_exothermic:
                funcs_to_save = [self.Tp, self.teq, self.Q_heat]
            else:
                funcs_to_save = [self.T]
            
            for func in funcs_to_save:
                filename = f"{self.output_dir}/{func.name}_passo_{time_step_idx:04d}.xdmf"
                with io.XDMFFile(self.comm, filename, "w") as xdmf:
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(func, current_time)

        except Exception as e:
            if self.rank == 0: print(f"     -> ⚠️ Erro ao salvar resultados: {e}")
            
    def _finalize(self):
        """Finaliza a simulação e imprime um resumo."""
        if self.rank == 0:
            print("\n" + "="*80)
            print("✅ SIMULAÇÃO CONCLUÍDA")
            print("="*80)
            print(f"   - Resultados salvos em: '{self.output_dir}'")
            print(f"   - Log de execução salvo em: '{self.logger.log_file.name}'")


def main():
    """Ponto de entrada principal para a execução da simulação."""
    import sys
    from pathlib import Path
    
    # Verificar argumentos de linha de comando
    if len(sys.argv) != 2:
        print("❌ USO: python barragem-Gemini-R0.py <nome_do_caso>")
        print("    Exemplo: python barragem-Gemini-R0.py barragem2")
        return
    
    # Nome do caso e diretório
    caso_arg = sys.argv[1]
    pasta_caso = Path(caso_arg)
    if not pasta_caso.is_dir():
        print(f"❌ Diretório do caso não encontrado: {pasta_caso}")
        return
    
    # Extrair o nome do caso (última parte do caminho)
    caso = pasta_caso.name
    
    # Montar caminhos dos arquivos
    config_file = pasta_caso / f"{caso}.yaml"
    json_file = pasta_caso / f"{caso}-xdmf.json"
    
    # Verificar se os arquivos existem
    if not config_file.exists():
        print(f"❌ Arquivo YAML não encontrado: {config_file}")
        return
    if not json_file.exists():
        print(f"❌ Arquivo JSON não encontrado: {json_file}")
        print(f"   Lembre-se de rodar 'debug_camadas.py {caso_arg}' para gerar o JSON.")
        return
    
    print(f"🔧 Iniciando simulação para caso: {caso}")
    print(f"   Configuração: {config_file}")
    print(f"   Análise: {json_file}")
    
    sim = SimulacaoBarragem(str(config_file), str(json_file))
    sim.run()

if __name__ == "__main__":
    main()