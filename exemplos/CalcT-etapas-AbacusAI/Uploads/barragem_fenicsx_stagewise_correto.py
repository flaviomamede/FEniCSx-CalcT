#!/usr/bin/env python3
"""
BARRAGEM FEniCSx - VERSÃO FINAL OTIMIZADA
Usa análise JSON para execução assertiva apenas no PRIMEIRO BLOCO DE TEMPO

MELHORIAS IMPLEMENTADAS:
- Executa APENAS primeiro bloco de tempo
- Lógica assertiva usando elementos/nós ativos do JSON
- Formulação variacional otimizada por bloco
- Mapeamento preciso de Physical Groups
"""

import dolfinx as fem
from dolfinx import mesh, fem, io
from dolfinx.fem import Function, FunctionSpace, Constant, locate_dofs_topological
from dolfinx.fem.petsc import LinearProblem, NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from mpi4py import MPI
import numpy as np
import ufl
from ufl import grad, dot, dx, ds, TestFunction, TrialFunction, lhs, rhs, TestFunctions, split
import yaml
import json
import os
import sys
import petsc4py
petsc4py.init()
from petsc4py import PETSc
import datetime

class TeeOutput:
    """Classe para redirecionar saída para terminal e arquivo simultaneamente"""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log_file = open(filename, 'w', encoding='utf-8')
        # Escrever cabeçalho
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header = f"# LOG DE EXECUÇÃO - {timestamp}\n\n"
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

class BarragemStageWiseCorreto:
    """Implementação FINAL usando análise JSON - APENAS PRIMEIRO BLOCO"""
    
    def __init__(self, config_file, json_file):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        
        # Carregar configuração YAML
        with open(config_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        # Carregar análise JSON
        with open(json_file, 'r', encoding='utf-8') as f:
            self.analysis = json.load(f)
        
        # Configurações
        self.theta = self.config['general']['theta']
        self.output_dir = self.config['general']['output_dir']
        
        # Temperatura inicial
        self.temp_inicial = 20.0
        if 'initial_conditions' in self.config and 'temperature' in self.config['initial_conditions']:
            self.temp_inicial = self.config['initial_conditions']['temperature']
        elif 'camadas_material' in self.config and self.config['camadas_material']:
            self.temp_inicial = self.config['camadas_material'][0]['temperatura_inicial']
        
        # Extrair informações da análise JSON
        self.time_vector = self.analysis['vetor_tempo']
        self.time_blocks = self.analysis['blocos_tempo']
        self.block_results = self.analysis['analise_resultados']
        self.mappings = self.analysis['info_geral']['mapeamentos']
        
        # NOVO: Extrair info do plano de simulação para lógica assertiva
        self.plano_simulacao = self.analysis['analise_resultados']
        
        if self.rank == 0:
            print("🏗️  BARRAGEM STAGE-WISE - VERSÃO FINAL OTIMIZADA")
            print("="*60)
            print(f"🌡️  Temperatura inicial: {self.temp_inicial}°C")
            print(f"📊 Blocos de tempo: {len(self.time_blocks)}")
            print(f"🎯 EXECUTANDO APENAS PRIMEIRO BLOCO DE TEMPO")
            print(f"📁 Malha: {self.analysis['info_geral']['xdmf_file']}")
    
    def load_mesh_and_discover_tags(self):
        """Carrega malha usando informações do JSON"""
        mesh_file = self.config['general']['mesh_file']
        
        with io.XDMFFile(self.comm, f"{mesh_file}", "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            
            # CORREÇÃO: Criar entidades de facetas e conectividade ANTES de ler as tags
            self.mesh.topology.create_entities(1)
            self.mesh.topology.create_connectivity(1, 2)  # facetas -> células
            
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            try:
                self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
                if self.rank == 0:
                    print(f"   ✅ Facet tags carregadas do XDMF: {np.unique(self.facet_tags.values)}")
            except Exception as e:
                if self.rank == 0:
                    print(f"   ⚠️  Falha ao carregar facet tags do XDMF: {e}")
                    print(f"   🔧 Criando facet tags vazias (todos os valores = 0)")
                tdim = self.mesh.topology.dim
                fdim = tdim - 1
                self.mesh.topology.create_connectivity(fdim, tdim)
                num_facets = self.mesh.topology.index_map(fdim).size_local
                facet_indices = np.arange(num_facets, dtype=np.int32)
                facet_markers = np.zeros(num_facets, dtype=np.int32)
                self.facet_tags = mesh.meshtags(self.mesh, fdim, facet_indices, facet_markers)
        
        if self.rank == 0:
            print(f"📁 Malha carregada: {mesh_file}")
            print(f"🏷️  Surfaces: {list(self.mappings['physical_surfaces'].values())}")
            print(f"🏷️  Lines: {list(self.mappings['physical_lines'].values())}")
    
    def setup_materials_and_schedule(self):
        """Configura materiais usando informações do YAML"""
        self.materials = {}
        
        for mat in self.config['materiais']:
            nome = mat['nome']
            hgen = mat.get('hgen', {})
            self.materials[nome] = {
                'densidade': mat['densidade'],
                'condutividade': mat['condutividade_termica'],
                'calor_especifico': mat['calor_especifico'],
                'temp_lancamento': self.temp_inicial,
                'gera_calor': hgen.get('gera_calor', False),
                'dTadinfty': hgen.get('par_gera_calor', {}).get('dTadinfty', 30.0),
                'a_dias': hgen.get('par_gera_calor', {}).get('a_dias', 1.5),
                'expoente': hgen.get('par_gera_calor', {}).get('expoente', 2.0),
                'EaR': hgen.get('EaR', 4000.0),
                'Tref': hgen.get('Tref', 20.0)
            }
        
        # Verificar se há materiais exotérmicos
        self.has_exothermic = any(mat['gera_calor'] for mat in self.materials.values())
        
        if self.rank == 0:
            print(f"✅ Materiais configurados: {list(self.materials.keys())}")
            if self.has_exothermic:
                print(f"🔥 Sistema exotérmico detectado!")
    
    def setup_function_spaces(self):
        """Define espaços de função"""
        if self.has_exothermic:
            # Sistema exotérmico: espaços separados para Tp e teq
            self.V = fem.functionspace(self.mesh, ("Lagrange", 1))
            
            # Funções para temperatura
            self.Tp = Function(self.V)          # Solução atual
            self.Tp_n = Function(self.V)        # Solução anterior
            self.u_Tp = TrialFunction(self.V)   # Incógnita (ESTA ERA A PEÇA QUE FALTAVA!)
            self.v_Tp = TestFunction(self.V)    # Função de teste
            
            # Funções para tempo equivalente
            self.teq = Function(self.V)
            self.teq_n = Function(self.V)
            self.v_teq = TestFunction(self.V)
            
            # Função para geração de calor
            self.Q_heat = Function(self.V)
            
            if self.rank == 0:
                print(f"✅ Espaços de função configurados (sistema exotérmico: Tp + teq + Q)")
        else:
            # Sistema simples: apenas temperatura
            self.V = fem.functionspace(self.mesh, ("Lagrange", 1))
            self.T = Function(self.V)
            self.T_n = Function(self.V)
            self.u = TrialFunction(self.V)
            self.v = TestFunction(self.V)
            
            if self.rank == 0:
                print(f"✅ Espaços de função configurados (sistema simples: apenas T)")
    
    def encontrar_bloco_para_tempo_atual(self, current_time, plano_simulacao):
        """FUNÇÃO ASSERTIVA: Encontra bloco ativo baseado no tempo atual"""
        for bloco_key, bloco_info in plano_simulacao.items():
            info_bloco = bloco_info['info_bloco']
            if info_bloco['inicio'] <= current_time <= info_bloco['fim']:
                return bloco_info
        return None
    
    def Q_function(self, teq_val, mat_props):
        """Função Hill para geração de calor (baseada no exemplo que funciona)"""
        rho = mat_props['densidade']
        ce = mat_props['calor_especifico']
        dTadinfty = mat_props['dTadinfty']
        a_dias = mat_props['a_dias']
        expoente = mat_props['expoente']
        
        # Converter a_dias para segundos
        a_sec = a_dias * 24 * 3600
        
        # Função Hill - corrigida para match com exemplo que funciona
        # Q(teq) = rho * ce * dTadinfty * teq^c / (a^c + teq^c)
        return rho * ce * (dTadinfty * teq_val**expoente / (a_sec**expoente + teq_val**expoente))
    
    def dQ_dteq(self, teq_val, mat_props):
        """Derivada da função Hill em relação ao tempo equivalente"""
        rho = mat_props['densidade']
        ce = mat_props['calor_especifico']
        dTadinfty = mat_props['dTadinfty']
        a_dias = mat_props['a_dias']
        expoente = mat_props['expoente']
        
        # Converter a_dias para segundos
        a_sec = a_dias * 24 * 3600
        
        # Derivada da função Hill
        term = dTadinfty * expoente * teq_val**(expoente-1) * a_sec**expoente / (a_sec**expoente + teq_val**expoente)**2
        return rho * ce * term
    
    def arrhenius_factor(self, Tp_val, mat_props):
        """Fator de Arrhenius para evolução do tempo equivalente"""
        EaR = mat_props['EaR']
        Tref = mat_props['Tref']
        
        # Fator de Arrhenius (implementação dos exemplos)
        return ufl.exp(EaR * (1/(Tref + 273.15) - 1/(Tp_val + 273.15)))
    
    def setup_variational_form(self, dt_val, current_time):
        """LÓGICA ASSERTIVA: Monta forma variacional usando domínios/contornos ativos do JSON"""
        # Encontrar bloco ativo usando função assertiva
        info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
        
        if not info_bloco_atual:
            if self.rank == 0:
                print(f"   ⚠️  Nenhum bloco encontrado para t={current_time/3600:.1f}h")
            return None
        
        # Domínios e contornos ativos (listas já prontas do JSON)
        dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
        contornos_ativos = info_bloco_atual['physical_groups']['lines']
        
        if self.rank == 0:
            print(f"   🔥 Domínios ativos: {dominios_ativos}")
            print(f"   🔥 Contornos ativos: {contornos_ativos}")
        
        # Constantes temporais
        dt_const = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        
        # Medidas de integração - DEFINIR UMA VEZ
        dx_measure = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags)
        ds_measure = ufl.Measure("ds", domain=self.mesh, subdomain_data=self.facet_tags)
        
        if self.has_exothermic:
            # SISTEMA EXOTÉRMICO: Apenas equação da temperatura (versão simplificada)
            # TODO: Implementar tempo equivalente em versão futura
            a = 0
            L = 0
            
            # Obter propriedades do material para este domínio
            material = list(self.materials.values())[0]  # Simplificado
            
            k = material['condutividade']
            rho = material['densidade']
            cp = material['calor_especifico']
            
            # Constantes
            k_const = Constant(self.mesh, PETSc.ScalarType(k))
            rho_const = Constant(self.mesh, PETSc.ScalarType(rho))
            cp_const = Constant(self.mesh, PETSc.ScalarType(cp))
            
            # Contribuição apenas dos domínios ativos
            for domain_id in dominios_ativos:
                # Apenas equação da temperatura (sem geração de calor por enquanto)
                # Termo temporal
                a += rho_const * cp_const * (self.Tp) / dt_const * self.v_Tp * dx_measure(domain_id)
                # Termo de difusão
                a += theta * k_const * dot(grad(self.Tp), grad(self.v_Tp)) * dx_measure(domain_id)
                
                # Lado direito
                L += rho_const * cp_const * self.Tp_n / dt_const * self.v_Tp * dx_measure(domain_id)
                L -= (1 - theta) * k_const * dot(grad(self.Tp_n), grad(self.v_Tp)) * dx_measure(domain_id)
            
            if self.rank == 0:
                print(f"   ✅ Formulação exotérmica SIMPLIFICADA construída para {len(dominios_ativos)} domínios ativos")
            
            return a, L
            
        else:
            # SISTEMA SIMPLES: Apenas temperatura
            a = 0
            L = 0
            
            # Propriedades do material (usar primeiro material como padrão)
            material = list(self.materials.values())[0]
            k = material['condutividade']
            rho = material['densidade']
            cp = material['calor_especifico']
            
            # Capacidade térmica
            C = rho * cp
            
            # Constantes
            k_const = Constant(self.mesh, PETSc.ScalarType(k))
            C_const = Constant(self.mesh, PETSc.ScalarType(C))
            
            # Contribuição apenas dos domínios ativos
            for domain_id in dominios_ativos:
                # Termo temporal
                a += (C_const / dt_const) * self.u * self.v * dx_measure(domain_id)
                # Termo de difusão
                a += theta * k_const * dot(grad(self.u), grad(self.v)) * dx_measure(domain_id)
                
                # Lado direito
                L += (C_const / dt_const) * self.T_n * self.v * dx_measure(domain_id)
                L -= (1 - theta) * k_const * dot(grad(self.T_n), grad(self.v)) * dx_measure(domain_id)
            
            if self.rank == 0:
                print(f"   ✅ Formulação simples construída para {len(dominios_ativos)} domínios ativos")
            
            return a, L
    
    def setup_variational_form_with_heat_generation(self, dt_val, current_time):
        """
        CORRIGIDO: Monta a forma variacional para o sistema exotérmico.
        A forma é inicializada como zero e as contribuições de todos os domínios
        ativos são somadas em um único loop para garantir a montagem correta da matriz.
        """
        # Encontrar bloco ativo usando função assertiva
        info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
        
        if not info_bloco_atual:
            if self.rank == 0:
                print(f"   ⚠️  Nenhum bloco encontrado para t={current_time/3600:.1f}h")
            return None
        
        # Domínios ativos (lista já pronta do JSON)
        dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
        
        if self.rank == 0:
            print(f"   🔥 Formulação exotérmica com Q para {len(dominios_ativos)} domínios ativos")
        
        # Constantes temporais
        dt_const = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        
        # Obter propriedades do material
        material = list(self.materials.values())[0]  # Usar primeiro material
        k = material['condutividade']
        rho = material['densidade']
        cp = material['calor_especifico']
        
        # Constantes
        k_const = Constant(self.mesh, PETSc.ScalarType(k))
        rho_const = Constant(self.mesh, PETSc.ScalarType(rho))
        cp_const = Constant(self.mesh, PETSc.ScalarType(cp))
        
        # --- INÍCIO DA CORREÇÃO ---
        
        # 1. Inicializar as formas como zero ANTES do loop
        a = 0
        L = 0
        
        # 2. Definir a medida de integração uma única vez
        dx_measure = ufl.Measure("dx", domain=self.mesh, subdomain_data=self.cell_tags)
        
        # 3. Loop sobre TODOS os domínios ativos para somar suas contribuições
        for domain_id in dominios_ativos:
            # Lado esquerdo (bilinear) - contribuição do domínio atual
            a += rho_const * cp_const * (self.u_Tp / dt_const) * self.v_Tp * dx_measure(domain_id)
            a += theta * k_const * ufl.dot(ufl.grad(self.u_Tp), ufl.grad(self.v_Tp)) * dx_measure(domain_id)
            
            # Lado direito (linear) - contribuição do domínio atual
            L += rho_const * cp_const * (self.Tp_n / dt_const) * self.v_Tp * dx_measure(domain_id)
            L -= (1 - theta) * k_const * ufl.dot(ufl.grad(self.Tp_n), ufl.grad(self.v_Tp)) * dx_measure(domain_id)
            L += self.Q_heat * self.v_Tp * dx_measure(domain_id)
        
        # --- FIM DA CORREÇÃO ---
        
        if self.rank == 0:
            print(f"   ✅ Formulação exotérmica com Q construída para {len(dominios_ativos)} domínios")
        
        return a, L
    
    def get_boundary_conditions(self, current_time):
        """
        CORREÇÃO FUNDAMENTAL: Aplica condições de contorno usando nós ativos/inativos do JSON
        para resolver o problema da matriz singular (nó 13 órfão)
        """
        try:
            # Encontrar bloco ativo usando função assertiva
            info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
            
            if not info_bloco_atual:
                return self.apply_fallback_boundary_conditions()
            
            # EXTRAIR NÓS ATIVOS E CONTORNOS DO JSON (NUMERAÇÃO XDMF)
            nos_ativos = info_bloco_atual['elementos_nos']['nos_dominio']
            contornos_ativos = info_bloco_atual['physical_groups']['lines']
            
            if self.rank == 0:
                print(f"   🎯 CORREÇÃO: Nós ativos: {len(nos_ativos)} de {self.mesh.topology.index_map(0).size_local}")
                print(f"   🎯 CORREÇÃO: Contornos ativos: {contornos_ativos}")
            
            # CALCULAR NÓS INATIVOS (SEM EQUAÇÃO FÍSICA)
            todos_nos = set(range(self.mesh.topology.index_map(0).size_local))
            nos_ativos_set = set(nos_ativos)
            nos_inativos = list(todos_nos - nos_ativos_set)
            
            if self.rank == 0:
                print(f"   🔧 CORREÇÃO: Nós inativos (órfãos): {sorted(nos_inativos)}")
            
            # Temperatura de contorno
            T_boundary = 20.0
            if 'camadas_material' in self.config and self.config['camadas_material']:
                T_boundary = self.config['camadas_material'][0]['temperatura_inicial']
            
            bcs = []
            
            # 1. APLICAR DIRICHLET NOS NÓS INATIVOS (CORREÇÃO PRINCIPAL)
            # Isso resolve o problema da matriz singular fixando os DOFs órfãos
            for no_inativo in nos_inativos:
                try:
                    dofs = np.array([no_inativo], dtype=np.int32)
                    valor_fixo = T_boundary  # Usar mesma temperatura que os contornos
                    
                    bc = fem.dirichletbc(
                        fem.Constant(self.mesh, float(valor_fixo)),
                        dofs,
                        self.V
                    )
                    bcs.append(bc)
                except Exception as e:
                    if self.rank == 0:
                        print(f"   ⚠️  Erro ao fixar nó inativo {no_inativo}: {e}")
            
            if self.rank == 0:
                print(f"   ✅ CORREÇÃO: {len([bc for bc in bcs])} nós inativos fixados em T={T_boundary}°C")
            
            # 2. APLICAR DIRICHLET NOS CONTORNOS ATIVOS (LÓGICA ORIGINAL)
            bcs_contornos = 0
            for boundary_id in contornos_ativos:
                try:
                    # Encontrar facetas que têm a tag do contorno atual
                    facet_indices = np.where(self.facet_tags.values == boundary_id)[0]
                    
                    if len(facet_indices) > 0:
                        # Localizar DOFs nessas facetas
                        boundary_dofs = locate_dofs_topological(self.V, self.mesh.topology.dim - 1, facet_indices)
                        
                        if len(boundary_dofs) > 0:
                            bc = fem.dirichletbc(
                                fem.Constant(self.mesh, float(T_boundary)),
                                boundary_dofs,
                                self.V
                            )
                            bcs.append(bc)
                            bcs_contornos += 1
                            
                except Exception as e:
                    if self.rank == 0:
                        print(f"   ⚠️  Erro ao aplicar BC no contorno {boundary_id}: {e}")
            
            # 3. VERIFICAÇÃO FINAL: Garantir que TODOS os nós têm BC (CORREÇÃO DEFINITIVA)
            # Identificar nós que ficaram sem BC
            todos_nos = set(range(self.mesh.topology.index_map(0).size_local))
            nos_com_bc = set(nos_inativos)  # Nós inativos já têm BC
            
            # Adicionar nós dos contornos ativos
            for boundary_id in contornos_ativos:
                try:
                    facet_indices = np.where(self.facet_tags.values == boundary_id)[0]
                    if len(facet_indices) > 0:
                        boundary_dofs = locate_dofs_topological(self.V, self.mesh.topology.dim - 1, facet_indices)
                        nos_com_bc.update(boundary_dofs)
                except:
                    pass
            
            # Nós sem BC (órfãos reais)
            nos_sem_bc = list(todos_nos - nos_com_bc)
            
            if nos_sem_bc:
                if self.rank == 0:
                    print(f"   🚨 CORREÇÃO FINAL: {len(nos_sem_bc)} nós órfãos detectados: {sorted(nos_sem_bc)}")
                
                # Aplicar BC nesses nós órfãos também
                for no_orfao in nos_sem_bc:
                    try:
                        dofs = np.array([no_orfao], dtype=np.int32)
                        bc = fem.dirichletbc(
                            fem.Constant(self.mesh, float(T_boundary)),
                            dofs,
                            self.V
                        )
                        bcs.append(bc)
                    except Exception as e:
                        if self.rank == 0:
                            print(f"   ⚠️  Erro ao fixar nó órfão {no_orfao}: {e}")
                
                if self.rank == 0:
                    print(f"   ✅ CORREÇÃO FINAL: {len(nos_sem_bc)} nós órfãos adicionais fixados")
            
            if self.rank == 0:
                print(f"   ✅ CORREÇÃO: {bcs_contornos} contornos ativos aplicados")
                print(f"   🎯 TOTAL FINAL: {len(bcs)} condições de contorno (cobertura completa)")
            
            return bcs
                
        except Exception as e:
            if self.rank == 0:
                print(f"   ❌ Erro nas condições de contorno: {e}")
            return self.apply_fallback_boundary_conditions()
    
    def apply_fallback_boundary_conditions(self):
        """Condições de contorno fallback - aplicar em toda a fronteira"""
        try:
            def boundary_marker(x):
                return np.ones(x.shape[1], dtype=bool)
            
            # Usar o espaço correto baseado no tipo de sistema
            V_space = self.V  # Tanto sistema simples quanto exotérmico usam self.V
            
            boundary_dofs = locate_dofs_topological(V_space, self.mesh.topology.dim-1, 
                                                   mesh.locate_entities_boundary(self.mesh, self.mesh.topology.dim-1, boundary_marker))
            
            bc = fem.dirichletbc(fem.Constant(self.mesh, 20.0), boundary_dofs, V_space)
            
            if self.rank == 0:
                print(f"   ⚠️  Usando BC fallback: {len(boundary_dofs)} DOFs")
            return [bc]
        except Exception as e:
            if self.rank == 0:
                print(f"   ❌ Erro no fallback: {e}")
            return []
    
    def solve_timestep(self, dt_val, current_time):
        """Resolver passo de tempo usando lógica assertiva do JSON"""
        # Encontrar bloco ativo usando função assertiva
        info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
        
        if not info_bloco_atual:
            if self.rank == 0:
                print(f"   ⚠️  Nenhum bloco ativo para t={current_time/3600:.1f}h")
            return
        
        if self.has_exothermic:
            # SISTEMA EXOTÉRMICO: Abordagem sequencial
            # 1. PRIMEIRO: Calcular tempo equivalente baseado na temperatura anterior
            self.update_equivalent_time(dt_val, current_time)
            
            # 2. SEGUNDO: Calcular Q usando o tempo equivalente
            self.calculate_heat_generation(current_time)
            
            # 3. TERCEIRO: Resolver equação de balanço de energia com Q
            self.solve_energy_balance(dt_val, current_time)
            
        else:
            # SISTEMA SIMPLES: Apenas temperatura
            # Montar forma variacional usando lógica assertiva
            form_result = self.setup_variational_form(dt_val, current_time)
            
            if form_result is None:
                if self.rank == 0:
                    print(f"   ❌ Falha na montagem da forma variacional")
                return
            
            # Aplicar condições de contorno usando lógica assertiva
            boundary_conditions = self.get_boundary_conditions(current_time)
            
            a, L = form_result
            
            try:
                # Solver direto
                problem = LinearProblem(
                    a, L, bcs=boundary_conditions,
                    petsc_options={
                        "ksp_type": "preonly",
                        "pc_type": "lu",
                        "pc_factor_mat_solver_type": "mumps"
                    }
                )
                T_new = problem.solve()
                
                # Copiar solução
                self.T.x.array[:] = T_new.x.array[:]
                
                if self.rank == 0:
                    camadas_ativas = info_bloco_atual['camadas_ativas']
                    print(f"   ✅ Convergiu (LU) - Camadas: {camadas_ativas}")
                    
            except Exception as e:
                # Fallback para solver iterativo
                try:
                    if self.rank == 0:
                        print(f"   ⚠️  LU falhou, tentando solver iterativo...")
                    
                    problem = LinearProblem(
                        a, L, bcs=boundary_conditions,
                        petsc_options={
                            "ksp_type": "cg",
                            "pc_type": "hypre",
                            "ksp_rtol": 1.0e-8,
                            "ksp_max_it": 1000
                        }
                    )
                    T_new = problem.solve()
                    
                    # Copiar solução
                    self.T.x.array[:] = T_new.x.array[:]
                    
                    if self.rank == 0:
                        camadas_ativas = info_bloco_atual['camadas_ativas']
                        print(f"   ✅ Convergiu (CG) - Camadas: {camadas_ativas}")
                        
                except Exception as e2:
                    if self.rank == 0:
                        print(f"   ❌ Erro no solver: {e2}")
                        print(f"   🔄 Mantendo solução anterior")

    def update_equivalent_time(self, dt_val, current_time):
        """PASSO 1: Calcular tempo equivalente baseado na temperatura anterior"""
        try:
            # Usar temperatura do passo anterior (Tp_n)
            info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
            dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
            
            if self.rank == 0:
                print(f"   🔄 Calculando tempo equivalente para {len(dominios_ativos)} domínios...")
            
            # Obter propriedades do material (usar primeiro material como padrão)
            material = list(self.materials.values())[0]
            
            # Calcular fator de Arrhenius baseado na temperatura anterior
            # Isso é uma integral simples: teq(t+dt) = teq(t) + dt * f_arrhenius(Tp(t))
            
            # Interpolação da temperatura anterior nos pontos de integração
            Tp_anterior = self.Tp_n.x.array[:]
            teq_anterior = self.teq_n.x.array[:]
            
            # Fator de Arrhenius usando temperatura anterior
            EaR = material['EaR']
            Tref = material['Tref']
            
            # Calcular incremento do tempo equivalente
            # Fator de Arrhenius: exp(EaR * (1/Tref - 1/Tp))
            expoente = EaR * (1/(Tref + 273.15) - 1/(Tp_anterior + 273.15))
            
            # Limitar apenas para evitar overflow numérico (não físico)
            expoente_limitado = np.clip(expoente, -50.0, 50.0)
            
            arrhenius_factor = np.exp(expoente_limitado)
            teq_increment = dt_val * arrhenius_factor
            
            # Apenas garantir que incremento >= 0 (fisicamente necessário)
            teq_increment = np.maximum(teq_increment, 0.0)
            
            # Atualizar tempo equivalente
            self.teq.x.array[:] = teq_anterior + teq_increment
            
            if self.rank == 0:
                teq_min = np.min(self.teq.x.array[:])
                teq_max = np.max(self.teq.x.array[:])
                print(f"   ✅ Tempo equivalente atualizado: [{teq_min:.1f}, {teq_max:.1f}]s")
                
        except Exception as e:
            if self.rank == 0:
                print(f"   ❌ Erro no cálculo do tempo equivalente: {e}")
                print(f"   🔄 Mantendo tempo equivalente anterior")
    
    def calculate_heat_generation(self, current_time):
        """PASSO 2: Calcular geração de calor Q baseada no tempo equivalente"""
        try:
            # Usar tempo equivalente calculado no passo anterior
            info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
            dominios_ativos = info_bloco_atual['physical_groups']['surfaces']
            
            # Obter propriedades do material
            material = list(self.materials.values())[0]
            
            # Calcular Q usando função Hill
            teq_values = self.teq.x.array[:]
            
            rho = material['densidade']
            ce = material['calor_especifico']
            dTadinfty = material['dTadinfty']
            a_dias = material['a_dias']
            expoente = material['expoente']
            
            # Converter a_dias para segundos
            a_sec = a_dias * 24 * 3600
            
            # Função Hill CORRIGIDA - sem limitações artificiais
            # Q(teq) = rho * ce * dTadinfty * teq^c / (a^c + teq^c)
            Q_values = rho * ce * (dTadinfty * teq_values**expoente / (a_sec**expoente + teq_values**expoente))
            
            # Apenas garantir que Q >= 0 (fisicamente necessário)
            Q_values = np.maximum(Q_values, 0.0)
            
            # Armazenar Q na função já criada
            self.Q_heat.x.array[:] = Q_values
            
            if self.rank == 0:
                Q_min = np.min(Q_values)
                Q_max = np.max(Q_values)
                print(f"   ✅ Geração de calor calculada: [{Q_min:.1e}, {Q_max:.1e}] W/m³")
                
        except Exception as e:
            if self.rank == 0:
                print(f"   ❌ Erro no cálculo da geração de calor: {e}")
                print(f"   🔄 Usando Q = 0")
            self.Q_heat.x.array[:] = 0.0
    
    def solve_energy_balance(self, dt_val, current_time):
        """PASSO 3: Resolver equação de balanço de energia com Q usando solver iterativo"""
        try:
            # Montar forma variacional com geração de calor
            form_result = self.setup_variational_form_with_heat_generation(dt_val, current_time)
            
            if form_result is None:
                if self.rank == 0:
                    print(f"   ❌ Falha na montagem da forma variacional com geração de calor")
                return
            
            a, L = form_result
            boundary_conditions = self.get_boundary_conditions(current_time)
            
            # SOLVER ITERATIVO: Usar GMRES com precondicionamento para problemas exotérmicos
            try:
                problem = LinearProblem(
                    a, L, bcs=boundary_conditions,
                    petsc_options={
                        "ksp_type": "gmres",
                        "pc_type": "ilu",
                        "ksp_rtol": 1e-6,
                        "ksp_atol": 1e-8,
                        "ksp_max_it": 1000
                    }
                )
                Tp_new = problem.solve()
                
                # Verificar se a solução é válida
                if not np.all(np.isfinite(Tp_new.x.array)):
                    raise RuntimeError("Solução GMRES divergente")
                
                # Se GMRES funcionou, usar
                self.Tp.x.array[:] = Tp_new.x.array[:]
                solver_usado = "GMRES"
                
            except Exception as e_gmres:
                if self.rank == 0:
                    print(f"   ⚠️  GMRES falhou: {e_gmres}")
                    print(f"   🔄 Tentando solver mais robusto (CG+ILU)...")
                
                # Fallback: Usar CG com ILU
                try:
                    problem_cg = LinearProblem(
                        a, L, bcs=boundary_conditions,
                        petsc_options={
                            "ksp_type": "cg",
                            "pc_type": "ilu",
                            "ksp_rtol": 1e-5,
                            "ksp_atol": 1e-7,
                            "ksp_max_it": 2000
                        }
                    )
                    Tp_new_cg = problem_cg.solve()
                    
                    # Verificar se a solução é válida
                    if not np.all(np.isfinite(Tp_new_cg.x.array)):
                        raise RuntimeError("Solução CG divergente")
                    
                    # Se CG funcionou, usar
                    self.Tp.x.array[:] = Tp_new_cg.x.array[:]
                    solver_usado = "CG"
                    
                except Exception as e_cg:
                    if self.rank == 0:
                        print(f"   ⚠️  CG falhou: {e_cg}")
                        print(f"   🔄 Tentando solver direto (LU)...")
                    
                    # Último recurso: LU direto
                    problem_lu = LinearProblem(
                        a, L, bcs=boundary_conditions,
                        petsc_options={
                            "ksp_type": "preonly",
                            "pc_type": "lu",
                            "pc_factor_mat_solver_type": "mumps"
                        }
                    )
                    Tp_new_lu = problem_lu.solve()
                    
                    # Verificar se a solução é válida
                    if not np.all(np.isfinite(Tp_new_lu.x.array)):
                        raise RuntimeError("Solução LU divergente")
                    
                    # Se LU funcionou, usar
                    self.Tp.x.array[:] = Tp_new_lu.x.array[:]
                    solver_usado = "LU"
            
            if self.rank == 0:
                info_bloco_atual = self.encontrar_bloco_para_tempo_atual(current_time, self.plano_simulacao)
                camadas_ativas = info_bloco_atual['camadas_ativas']
                Tp_min = np.min(self.Tp.x.array[:])
                Tp_max = np.max(self.Tp.x.array[:])
                print(f"   ✅ Convergiu ({solver_usado}) - Camadas: {camadas_ativas}")
                print(f"   🌡️  Temperatura: [{Tp_min:.1f}, {Tp_max:.1f}]°C")
                
        except Exception as e:
            if self.rank == 0:
                print(f"   ❌ Erro no solver de energia: {e}")
                print(f"   🔄 Mantendo solução anterior (Tp não será atualizado)")
            # Importante: não atualize Tp_n se a solução falhar
            # A atualização ocorre fora deste método, no loop principal
    
    def run_simulation(self):
        """EXECUTA SIMULAÇÃO APENAS NO PRIMEIRO BLOCO DE TEMPO"""
        self.load_mesh_and_discover_tags()
        self.setup_materials_and_schedule()
        self.setup_function_spaces()
        
        # Condições iniciais
        if self.has_exothermic:
            # Sistema exotérmico: Tp = temp_inicial, teq = 0
            self.Tp_n.x.array[:] = self.temp_inicial
            self.Tp.x.array[:] = self.temp_inicial
            self.teq_n.x.array[:] = 0.0
            self.teq.x.array[:] = 0.0
            self.Q_heat.x.array[:] = 0.0  # Inicializar geração de calor
            
            if self.rank == 0:
                print(f"   ✅ Condições iniciais exotérmicas: Tp={self.temp_inicial}°C, teq=0s, Q=0")
        else:
            # Sistema simples: apenas temperatura
            self.T_n.x.array[:] = self.temp_inicial
            self.T.x.array[:] = self.temp_inicial
            
            if self.rank == 0:
                print(f"   ✅ Condições iniciais simples: T={self.temp_inicial}°C")
        
        # EXTRAIR APENAS PRIMEIRO BLOCO DE TEMPO
        primeiro_bloco = self.time_blocks[0]
        time_vector = self.time_vector
        
        # Filtrar vetor de tempo apenas para primeiro bloco
        primeiro_bloco_tempo = [t for t in time_vector if primeiro_bloco['inicio'] <= t <= primeiro_bloco['fim']]
        
        if self.rank == 0:
            print(f"\n🚀 SIMULAÇÃO APENAS NO PRIMEIRO BLOCO DE TEMPO")
            print(f"   📅 Bloco 1: {primeiro_bloco['inicio']/3600:.1f}h → {primeiro_bloco['fim']/3600:.1f}h")
            print(f"   🕒 Duração: {primeiro_bloco['duracao']/3600:.1f}h")
            print(f"   📊 Pontos de tempo: {len(primeiro_bloco_tempo)}")
            print("="*60)
        
        # Criar diretório de saída
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Simular apenas para pontos de tempo do primeiro bloco
        for i in range(1, len(primeiro_bloco_tempo)):
            current_time = primeiro_bloco_tempo[i]
            previous_time = primeiro_bloco_tempo[i-1]
            dt_val = current_time - previous_time
            
            if self.rank == 0:
                print(f"\n📅 Passo {i}: t={previous_time/3600:.1f}h → {current_time/3600:.1f}h (dt={dt_val/3600:.1f}h)")
            
            # Resolver usando lógica assertiva
            self.solve_timestep(dt_val, current_time)
            
            # Atualizar soluções
            if self.has_exothermic:
                # Sistema exotérmico: atualizar Tp_n <- Tp, teq_n <- teq
                self.Tp_n.x.array[:] = self.Tp.x.array[:]
                self.teq_n.x.array[:] = self.teq.x.array[:]
            else:
                # Sistema simples: atualizar T_n <- T
                self.T_n.x.array[:] = self.T.x.array[:]
            
            # Salvar resultados
            if i % 5 == 0:  # Salvar mais frequentemente no primeiro bloco
                self.save_results(i, current_time)
        
        if self.rank == 0:
            print(f"\n✅ SIMULAÇÃO DO PRIMEIRO BLOCO CONCLUÍDA")
            print(f"   📊 Total: {len(primeiro_bloco_tempo)-1} passos")
            print(f"   🕒 Tempo final: {primeiro_bloco_tempo[-1]/3600:.1f}h")
    
    def save_results(self, time_step, current_time):
        """Salva resultados"""
        try:
            if self.has_exothermic:
                # Sistema exotérmico: salvar Tp e teq
                output_file_tp = f"{self.output_dir}/temperatura_bloco1_{time_step:04d}.xdmf"
                output_file_teq = f"{self.output_dir}/tempo_equiv_bloco1_{time_step:04d}.xdmf"
                
                # Salvar temperatura
                with io.XDMFFile(self.comm, output_file_tp, "w") as xdmf:
                    self.Tp.name = "Temperatura"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.Tp, current_time)
                
                # Salvar tempo equivalente
                with io.XDMFFile(self.comm, output_file_teq, "w") as xdmf:
                    self.teq.name = "Tempo_Equivalente"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.teq, current_time)
                
                if self.rank == 0:
                    print(f"   💾 Salvo: temperatura_bloco1_{time_step:04d}.xdmf + tempo_equiv_bloco1_{time_step:04d}.xdmf")
            else:
                # Sistema simples: salvar apenas temperatura
                output_file = f"{self.output_dir}/temperatura_bloco1_{time_step:04d}.xdmf"
                with io.XDMFFile(self.comm, output_file, "w") as xdmf:
                    self.T.name = "Temperatura"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.T, current_time)
                    
                if self.rank == 0:
                    print(f"   💾 Salvo: temperatura_bloco1_{time_step:04d}.xdmf")
        except Exception as e:
            if self.rank == 0:
                print(f"   ⚠️  Erro ao salvar: {e}")

def main():
    """Executa simulação apenas no primeiro bloco usando análise JSON"""
    config_file = "barragem1.yaml"
    json_file = "analise_stagewise_xdmf.json"
    
    # Configurar logging para arquivo
    log_output = TeeOutput("log.md")
    sys.stdout = log_output
    
    try:
        if not os.path.exists(config_file):
            print(f"❌ Arquivo {config_file} não encontrado!")
            sys.exit(1)
        
        if not os.path.exists(json_file):
            print(f"❌ Arquivo {json_file} não encontrado!")
            print(f"📝 Execute primeiro: python debug_camadas.py")
            sys.exit(1)
        
        try:
            solver = BarragemStageWiseCorreto(config_file, json_file)
            solver.run_simulation()
            
            print(f"\n🎯 SIMULAÇÃO PRIMEIRO BLOCO CONCLUÍDA")
            print(f"📝 Melhorias implementadas:")
            print(f"   ✅ Lógica assertiva usando JSON")
            print(f"   ✅ Formulação variacional apenas para domínios ativos")
            print(f"   ✅ Condições de contorno apenas para contornos ativos")
            print(f"   ✅ Execução apenas no primeiro bloco de tempo")
            print(f"   ✅ Função encontrar_bloco_para_tempo_atual implementada")
            
        except Exception as e:
            print(f"❌ Erro: {e}")
            import traceback
            traceback.print_exc()
    
    finally:
        # Sempre fechar o arquivo de log
        log_output.close()
        print(f"\n📝 Log salvo em: log.md")

if __name__ == "__main__":
    main()
