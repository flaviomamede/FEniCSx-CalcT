#!/usr/bin/env python3
"""
Simulação Térmica FEniCSx - GENÉRICO
Lê Physical Groups dinamicamente e implementa etapas construtivas
"""

import numpy as np
import yaml
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
from dolfinx import mesh, fem, io
from dolfinx.fem import FunctionSpace, Function, Constant
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
import ufl
from ufl import grad, dot, dx, ds, inner, TestFunction, TestFunctions, lhs, rhs
import os
import h5py

class GenericThermalFEniCSx:
    def __init__(self, yaml_file, xdmf_file):
        """Inicializa solver genérico"""
        self.yaml_file = yaml_file
        self.xdmf_file = xdmf_file
        self.h5_file = xdmf_file.replace('.xdmf', '.h5')
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        
        # Carregar configuração
        self.load_configuration()
        
        # Carregar malha e descobrir Physical Groups
        self.load_mesh_and_discover_tags()
        
        # Mapear dinamicamente Physical Groups -> YAML
        self.map_physical_groups_to_yaml()
        
        # Configurar problema (materials primeiro para detectar exotérmico)
        self.setup_materials_and_layers()
        self.setup_function_spaces()
        self.setup_initial_conditions()
        
    def load_configuration(self):
        """Carrega YAML"""
        with open(self.yaml_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self.tempo_final = self.config['general']['tempo_final']
        self.delta_t = self.config['general']['delta_t']
        self.delta_t_refinado = self.config['general']['delta_t_refinado']
        self.theta = self.config['general']['theta']
        self.output_dir = self.config['general']['output_dir']
        
        if self.rank == 0:
            os.makedirs(self.output_dir, exist_ok=True)
            print(f"📋 Configuração carregada: {self.tempo_final/86400:.1f} dias")
    
    def load_mesh_and_discover_tags(self):
        """Carrega malha e descobre automaticamente os Physical Groups"""
        # Carregar malha
        with io.XDMFFile(self.comm, self.xdmf_file, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            
            # Criar entidades de facetas antes de ler tags
            self.mesh.topology.create_connectivity(self.mesh.topology.dim-1, 0)
            
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        
        # Descobrir Physical Groups automaticamente
        self.discovered_cell_tags = np.unique(self.cell_tags.values)
        self.discovered_facet_tags = np.unique(self.facet_tags.values)
        
        if self.rank == 0:
            print(f"✅ Malha carregada: {self.mesh.topology.index_map(0).size_global} nós")
            print(f"🔍 Physical Groups descobertos:")
            print(f"   Células (volumes): {self.discovered_cell_tags}")
            print(f"   Facetas (contornos): {self.discovered_facet_tags}")
        
        # Ler nomes dos Physical Groups do arquivo H5 (se disponível)
        self.read_physical_group_names()
    
    def read_physical_group_names(self):
        """Tenta ler nomes dos Physical Groups do arquivo H5"""
        self.cell_tag_names = {}
        self.facet_tag_names = {}
        
        try:
            # Esta é uma implementação básica - pode precisar ser ajustada
            # dependendo de como o gmshio salva os nomes
            if self.rank == 0:
                print("🔍 Tentando ler nomes dos Physical Groups...")
                
                # Por enquanto, criar nomes genéricos baseados nos IDs
                for tag in self.discovered_cell_tags:
                    self.cell_tag_names[tag] = f"volume_{tag}"
                
                for tag in self.discovered_facet_tags:
                    self.facet_tag_names[tag] = f"boundary_{tag}"
                    
                print(f"📝 Nomes das células: {self.cell_tag_names}")
                print(f"📝 Nomes das facetas: {self.facet_tag_names}")
                
        except Exception as e:
            if self.rank == 0:
                print(f"⚠️  Não foi possível ler nomes: {e}")
    
    def map_physical_groups_to_yaml(self):
        """Mapeia automaticamente Physical Groups descobertos para configuração YAML"""
        # Mapear camadas de material (cell_tags) para materiais
        self.material_mapping = {}
        self.active_layers = {}  # Controla quais camadas estão ativas
        
        for camada_mat in self.config['camadas_material']:
            nome_camada = camada_mat['nome']
            material = camada_mat['material']
            
            # Extrair ID da camada
            if 'camada_material_' in nome_camada:
                camada_id = int(nome_camada.split('_')[-1])
                
                # Verificar se esta tag existe na malha
                if camada_id in self.discovered_cell_tags:
                    self.material_mapping[camada_id] = material
                    # Inicialmente, todas as camadas estão inativas
                    self.active_layers[camada_id] = False
        
        # Mapear contornos (facet_tags) para condições de contorno
        self.boundary_conditions = {}
        
        for contorno in self.config['contornos']:
            nome = contorno['nome']
            
            # Tentar encontrar a tag correspondente
            # Aqui você pode implementar lógica mais sofisticada de matching
            found_tag = None
            
            # Procurar por correspondência nos nomes descobertos
            for tag in self.discovered_facet_tags:
                tag_name = self.facet_tag_names.get(tag, f"boundary_{tag}")
                
                # Implementar matching inteligente aqui
                # Por enquanto, usar mapeamento direto baseado no conhecimento do .geo
                mapping_dict = {
                    'ISOLAMENTO_PERFEITO': 11,
                    'FUNDACAO_TOPO': 12,
                    'FACE_MONTANTE_1': 13,
                    'FACE_JUSANTE_1': 14,
                    'FACE_MONTANTE_2': 15,
                    'FACE_JUSANTE_2': 16,
                    'FACE_MONTANTE_3': 17,
                    'FACE_JUSANTE_3': 18,
                    'FACE_TOPO': 19,
                    'interface_1_2': 20,
                    'interface_2_3': 21
                }
                
                if nome in mapping_dict and tag == mapping_dict[nome]:
                    found_tag = tag
                    break
            
            if found_tag is not None:
                self.boundary_conditions[found_tag] = contorno
        
        if self.rank == 0:
            print(f"🗺️  Mapeamento materiais: {self.material_mapping}")
            print(f"🗺️  Condições de contorno mapeadas: {len(self.boundary_conditions)} encontradas")
    
    def setup_function_spaces(self):
        """Define espaços de função - detecta automaticamente se precisa sistema exotérmico"""
        # TEMPORÁRIO: Forçar sistema simples para debug
        # self.has_exothermic = any(
        #     mat.get('gera_calor', False) for mat in self.materials.values()
        # )
        self.has_exothermic = False  # TESTE
        
        if self.has_exothermic:
            # Sistema exotérmico: Temperatura + Tempo Equivalente (sistema misto simples)
            self.V = fem.functionspace(self.mesh, ("Lagrange", 1))
            
            # Criar duas funções separadas para Tp e teq (mais simples que elemento misto)
            self.T = Function(self.V)    # Temperatura
            self.T_n = Function(self.V)  # Temperatura anterior
            self.teq = Function(self.V)  # Tempo equivalente
            self.teq_n = Function(self.V) # Tempo equivalente anterior
            
            # Funções de teste
            self.v_T = TestFunction(self.V)   # Para temperatura
            self.v_teq = TestFunction(self.V) # Para tempo equivalente
            
            if self.rank == 0:
                print(f"✅ Sistema EXOTÉRMICO configurado (Tp + teq)")
        else:
            # Sistema simples: apenas Temperatura
            self.V = fem.functionspace(self.mesh, ("Lagrange", 1))
            self.T = Function(self.V)
            self.Tn = Function(self.V)
            self.v = TestFunction(self.V)
            
            if self.rank == 0:
                print(f"✅ Sistema SIMPLES configurado (apenas Temperatura)")
    
    def setup_materials_and_layers(self):
        """Configura materiais e camadas com birth/death"""
        # Processar materiais
        self.materials = {}
        for mat in self.config['materiais']:
            nome = mat['nome']
            self.materials[nome] = {
                'densidade': mat['densidade'],
                'condutividade': mat['condutividade_termica'],
                'calor_especifico': mat['calor_especifico'],
                'gera_calor': mat.get('hgen', {}).get('gera_calor', False),
                'termoactivation': mat.get('hgen', {}).get('termoactivation', False)
            }
            
            if 'hgen' in mat and 'par_gera_calor' in mat['hgen']:
                par = mat['hgen']['par_gera_calor']
                self.materials[nome]['Tad_inf'] = par['dTadinfty']
                self.materials[nome]['a_dias'] = par['a_dias']
                self.materials[nome]['expoente'] = par['expoente']
        
        # Processar camadas com birth/death
        self.layer_schedule = {}
        for camada in self.config['camadas']:
            nome = camada['nome']
            birth_time = camada['birth']
            death_time = camada.get('death', None)
            
            self.layer_schedule[nome] = {
                'birth': birth_time,
                'death': death_time
            }
        
        if self.rank == 0:
            print(f"✅ Materiais configurados: {list(self.materials.keys())}")
            print(f"⏰ Cronograma de camadas: {self.layer_schedule}")
    
    def setup_initial_conditions(self):
        """Define condições iniciais"""
        # Configurar parâmetros exotérmicos se necessário
        self.setup_exothermic_functions()
        
        # Obter temperatura inicial da configuração
        temp_default = 20.0
        if 'initial_conditions' in self.config and 'temperature' in self.config['initial_conditions']:
            temp_default = self.config['initial_conditions']['temperature']
        else:
            # Usar temperatura inicial da primeira camada material
            for camada_mat in self.config['camadas_material']:
                temp_default = camada_mat['temperatura_inicial']
                break
        
        if self.has_exothermic:
            # Sistema exotérmico: inicializar Tp e teq separadamente
            def initial_temp(x):
                return np.full(x.shape[1], temp_default)
            
            def initial_teq(x):
                return np.zeros(x.shape[1])  # Tempo equivalente inicial = 0
            
            # Inicializar temperatura
            self.T_n.interpolate(initial_temp)
            self.T.x.array[:] = self.T_n.x.array[:]
            
            # Inicializar tempo equivalente  
            self.teq_n.interpolate(initial_teq)
            self.teq.x.array[:] = self.teq_n.x.array[:]
            
            if self.rank == 0:
                print(f"✅ Condições iniciais EXOTÉRMICAS: Tp={temp_default}°C, teq=0s")
        else:
            # Sistema simples: apenas temperatura
            def initial_temp(x):
                return np.full(x.shape[1], temp_default)
            
            self.Tn.interpolate(initial_temp)
            self.T.x.array[:] = self.Tn.x.array[:]
            
            if self.rank == 0:
                print(f"✅ Condições iniciais SIMPLES: T={temp_default}°C")
    
    def setup_exothermic_functions(self):
        """Configura funções para geração de calor exotérmica"""
        if not self.has_exothermic:
            return
            
        # Parâmetros do problema exotérmico (do arquivo .pde)
        self.Tad_inf1 = 16.40021368590
        self.a1 = 0.42998953984 * 24 * 3600  # Convertendo para segundos
        self.c1 = 1.93329591299
        self.Tad_inf2 = 13.59978631410
        self.a2 = 1.34507960663 * 24 * 3600  # Convertendo para segundos
        self.c2 = 1.57098129317
        self.EaR = 4000  # Energia de ativação/R (K)
        
        if self.rank == 0:
            print(f"✅ Parâmetros exotérmicos configurados")
    
    def Q_function(self, teq_val, rho, ce):
        """Calcula geração de calor Q baseado no tempo equivalente"""
        return rho * ce * (
            self.Tad_inf1 * teq_val**self.c1 / (self.a1**self.c1 + teq_val**self.c1) + 
            self.Tad_inf2 * teq_val**self.c2 / (self.a2**self.c2 + teq_val**self.c2)
        )
    
    def dQ_dteq(self, teq_val, rho, ce):
        """Derivada de Q em relação ao tempo equivalente"""
        term1 = (self.Tad_inf1 * self.c1 * teq_val**(self.c1-1) * self.a1**self.c1 / 
                (self.a1**self.c1 + teq_val**self.c1)**2)
        term2 = (self.Tad_inf2 * self.c2 * teq_val**(self.c2-1) * self.a2**self.c2 / 
                (self.a2**self.c2 + teq_val**self.c2)**2)
        return rho * ce * (term1 + term2)
    
    def arrhenius_factor(self, Tp_val):
        """Fator de Arrhenius para evolução do tempo equivalente"""
        return ufl.exp(self.EaR * (1/298.15 - 1/(Tp_val + 273.15)))
    
    def update_active_layers(self, current_time):
        """Atualiza quais camadas estão ativas baseado no birth/death"""
        # Resetar estado das camadas
        for layer_id in self.active_layers:
            self.active_layers[layer_id] = False
        
        # Ativar camadas baseado no cronograma
        for layer_name, schedule in self.layer_schedule.items():
            birth_time = schedule['birth']
            death_time = schedule['death']
            
            # Verificar se camada deve estar ativa
            is_active = current_time >= birth_time
            if death_time is not None:
                is_active = is_active and current_time < death_time
            
            # Encontrar camadas de material correspondentes
            for camada_mat in self.config['camadas_material']:
                if camada_mat['camada'] == layer_name:
                    nome_camada = camada_mat['nome']
                    if 'camada_material_' in nome_camada:
                        camada_id = int(nome_camada.split('_')[-1])
                        if camada_id in self.active_layers:
                            self.active_layers[camada_id] = is_active
        
        # Log das camadas ativas
        if self.rank == 0:
            active_list = [lid for lid, active in self.active_layers.items() if active]
            print(f"⚡ Camadas ativas no tempo {current_time/3600:.1f}h: {active_list}")
    
    def get_material_properties(self, domain_id):
        """Obtém propriedades do material para domínio específico"""
        if domain_id in self.material_mapping:
            material_name = self.material_mapping[domain_id]
            return self.materials[material_name]
        else:
            # Material padrão
            return list(self.materials.values())[0]
    
    def get_active_boundaries(self, current_time):
        """Determina contornos ativos baseado nos domínios ativos e regras construtivas"""
        # Obter domínios ativos atuais
        dominios_ativos = [lid for lid, active in self.active_layers.items() if active]
        contornos_ativos = []
        
        # REGRAS DE ATIVAÇÃO (baseadas no debug_contornos.py)
        contornos_config = {
            11: {'conecta': 'fundacao'},      # ISOLAMENTO_PERFEITO  
            12: {'conecta': 'fundacao'},      # FUNDACAO_TOPO
            13: {'conecta': 'camada_1'},      # FACE_MONTANTE_1
            14: {'conecta': 'camada_1'},      # FACE_JUSANTE_1  
            15: {'conecta': 'camada_2'},      # FACE_MONTANTE_2
            16: {'conecta': 'camada_2'},      # FACE_JUSANTE_2
            17: {'conecta': 'camada_3'},      # FACE_MONTANTE_3
            18: {'conecta': 'camada_3'},      # FACE_JUSANTE_3
            19: {'conecta': 'topo_ativo'},    # FACE_TOPO (muda conforme última camada)
            20: {'conecta': 'interface_1_2'}, # interface 1-2 (morre quando camada 2 nasce)
            21: {'conecta': 'interface_2_3'}  # interface 2-3 (morre quando camada 3 nasce)
        }
        
        for tag, info in contornos_config.items():
            if tag not in self.boundary_conditions:
                continue  # Contorno não configurado no YAML
                
            conecta = info['conecta']
            ativo = False
            
            if conecta == 'fundacao':
                # Contornos da fundação: sempre ativos se fundação ativa
                ativo = any(d in [1,2,3,4] for d in dominios_ativos)
                
            elif conecta == 'camada_1':
                # Contornos externos da camada 1: ativos se camada 1 ativa
                ativo = any(d in [5,6] for d in dominios_ativos)
                
            elif conecta == 'camada_2':
                # Contornos externos da camada 2: ativos se camada 2 ativa
                ativo = any(d in [7,8] for d in dominios_ativos)
                
            elif conecta == 'camada_3':
                # Contornos externos da camada 3: ativos se camada 3 ativa
                ativo = any(d in [9,10] for d in dominios_ativos)
                
            elif conecta == 'topo_ativo':
                # Topo: muda conforme última camada ativa
                if any(d in [9,10] for d in dominios_ativos):
                    ativo = True  # topo da camada 3
                elif any(d in [7,8] for d in dominios_ativos):
                    ativo = True  # topo da camada 2  
                elif any(d in [5,6] for d in dominios_ativos):
                    ativo = True  # topo da camada 1
                    
            elif conecta == 'interface_1_2':
                # Interface 1-2: ATIVA só quando camada 1 existe mas camada 2 NÃO
                tem_camada1 = any(d in [5,6] for d in dominios_ativos)
                tem_camada2 = any(d in [7,8] for d in dominios_ativos)
                ativo = tem_camada1 and not tem_camada2
                
            elif conecta == 'interface_2_3':
                # Interface 2-3: ATIVA só quando camada 2 existe mas camada 3 NÃO
                tem_camada2 = any(d in [7,8] for d in dominios_ativos)
                tem_camada3 = any(d in [9,10] for d in dominios_ativos)
                ativo = tem_camada2 and not tem_camada3
            
            if ativo:
                contornos_ativos.append(tag)
        
        if self.rank == 0:
            print(f"🔗 Contornos ativos no tempo {current_time/3600:.1f}h: {contornos_ativos}")
        
        return contornos_ativos
    
    def setup_variational_form(self, dt_val, current_time):
        """Formulação variacional - simples ou exotérmica"""
        dt = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        
        # Medidas de integração
        dx_tags = dx(domain=self.mesh, subdomain_data=self.cell_tags)
        ds_tags = ds(domain=self.mesh, subdomain_data=self.facet_tags)
        
        if self.has_exothermic:
            return self.setup_exothermic_form(dt, theta, dx_tags, ds_tags, current_time)
        else:
            return self.setup_simple_form(dt, theta, dx_tags, ds_tags, current_time)
    
    def setup_simple_form(self, dt, theta, dx_tags, ds_tags, current_time):
        """Formulação simples (apenas temperatura)"""
        F = 0
        
        # Processar apenas domínios ativos
        for domain_id in self.discovered_cell_tags:
            # Verificar se camada está ativa
            if domain_id not in self.active_layers or not self.active_layers[domain_id]:
                continue
            
            # Obter propriedades do material
            mat_props = self.get_material_properties(domain_id)
            
            # Constantes do material
            rho = Constant(self.mesh, PETSc.ScalarType(mat_props['densidade']))
            k = Constant(self.mesh, PETSc.ScalarType(mat_props['condutividade']))
            cp = Constant(self.mesh, PETSc.ScalarType(mat_props['calor_especifico']))
            
            # Formulação de Crank-Nicolson
            F += rho * cp * (self.T - self.Tn) / dt * self.v * dx_tags(domain_id)
            
            T_theta = theta * self.T + (1 - theta) * self.Tn
            F += k * dot(grad(T_theta), grad(self.v)) * dx_tags(domain_id)
        
        # Aplicar APENAS condições de contorno dos contornos ATIVOS 
        contornos_ativos = self.get_active_boundaries(current_time)
        
        for boundary_tag in contornos_ativos:
            if boundary_tag in self.boundary_conditions:
                bc_config = self.boundary_conditions[boundary_tag]
                if bc_config['tipo'] == 'conveccao':
                    h_val = bc_config['h']
                    T_ext_val = bc_config['t_ext']
                    
                    h = Constant(self.mesh, PETSc.ScalarType(h_val))
                    T_ext = Constant(self.mesh, PETSc.ScalarType(T_ext_val))
                    
                    T_boundary = theta * self.T + (1 - theta) * self.Tn
                    F += h * (T_boundary - T_ext) * self.v * ds_tags(boundary_tag)
        
        return F
    
    def setup_exothermic_form(self, dt, theta, dx_tags, ds_tags):
        """Formulação exotérmica (Tp + teq acoplados)"""
        F_Tp = 0
        F_teq = 0
        
        # Processar apenas domínios ativos
        for domain_id in self.discovered_cell_tags:
            # Verificar se camada está ativa
            if domain_id not in self.active_layers or not self.active_layers[domain_id]:
                continue
            
            # Obter propriedades do material
            mat_props = self.get_material_properties(domain_id)
            
            # Constantes do material
            rho_val = mat_props['densidade']
            k_val = mat_props['condutividade']
            cp_val = mat_props['calor_especifico']
            
            rho = Constant(self.mesh, PETSc.ScalarType(rho_val))
            k = Constant(self.mesh, PETSc.ScalarType(k_val))
            cp = Constant(self.mesh, PETSc.ScalarType(cp_val))
            
            if mat_props['gera_calor']:
                # ===== EQUAÇÃO DE TEMPERATURA =====
                # Termo temporal
                F_Tp += rho * cp * (self.T - self.T_n) / dt * self.v_T * dx_tags(domain_id)
                
                # Termo de difusão
                T_theta = theta * self.T + (1 - theta) * self.T_n
                F_Tp += k * dot(grad(T_theta), grad(self.v_T)) * dx_tags(domain_id)
                
                # Termo exotérmico (derivada de Q)
                dQ_val = self.dQ_dteq(self.teq_n, rho_val, cp_val)
                F_Tp -= dQ_val * (self.teq - self.teq_n) / dt * self.v_T * dx_tags(domain_id)
                
                # ===== EQUAÇÃO DO TEMPO EQUIVALENTE =====
                # Evolução do tempo equivalente
                F_teq += (self.teq - self.teq_n) / dt * self.v_teq * dx_tags(domain_id)
                
                # Termo de Arrhenius
                arrhenius = self.arrhenius_factor(self.T_n)
                F_teq -= arrhenius * self.v_teq * dx_tags(domain_id)
            else:
                # Material sem geração de calor - apenas temperatura
                F_Tp += rho * cp * (self.T - self.T_n) / dt * self.v_T * dx_tags(domain_id)
                T_theta = theta * self.T + (1 - theta) * self.T_n
                F_Tp += k * dot(grad(T_theta), grad(self.v_T)) * dx_tags(domain_id)
                
                # Tempo equivalente permanece constante
                F_teq += (self.teq - self.teq_n) * self.v_teq * dx_tags(domain_id)
        
        # Aplicar condições de contorno para temperatura
        for boundary_tag, bc_config in self.boundary_conditions.items():
            if bc_config['tipo'] == 'conveccao':
                h_val = bc_config['h']
                T_ext_val = bc_config['t_ext']
                
                h = Constant(self.mesh, PETSc.ScalarType(h_val))
                T_ext = Constant(self.mesh, PETSc.ScalarType(T_ext_val))
                
                T_boundary = theta * self.T + (1 - theta) * self.T_n
                F_Tp += h * (T_boundary - T_ext) * self.v_T * ds_tags(boundary_tag)
        
        return F_Tp + F_teq
    
    def solve_timestep(self, dt_val, current_time):
        """Resolve passo de tempo"""
        # Atualizar camadas ativas
        self.update_active_layers(current_time)
        
        # Verificar se há camadas ativas
        if not any(self.active_layers.values()):
            if self.rank == 0:
                print(f"⚠️  Nenhuma camada ativa no tempo {current_time/3600:.1f}h")
            return
        
        bcs = []  # Condições Dirichlet (se houver)
        
        if self.has_exothermic:
            # Sistema exotérmico: resolver Tp e teq sequencialmente
            self.solve_exothermic_timestep(dt_val, current_time, bcs)
        else:
            # Sistema simples: resolver apenas temperatura
            F = self.setup_variational_form(dt_val, current_time)
            problem = NonlinearProblem(F, self.T, bcs)
            solver = NewtonSolver(self.comm, problem)
            
            # Configuração mais robusta
            solver.convergence_criterion = "residual"
            solver.rtol = 1e-3  # Tolerância mais relaxada
            solver.atol = 1e-6  # Tolerância absoluta
            solver.max_it = 20  # Menos iterações
            
            try:
                n_iterations, converged = solver.solve(self.T)
                if not converged:
                    # Tentar com método linear direto para sistemas mais simples
                    if self.rank == 0:
                        print(f"⚠️  Newton não convergiu, tentando solver linear...")
                    
                    # Reset para valor anterior e usar um passo menor
                    from dolfinx.fem.petsc import LinearProblem
                    try:
                        linear_problem = LinearProblem(lhs(F), rhs(F), self.T, bcs)
                        linear_problem.solve()
                        if self.rank == 0:
                            print(f"✅ Solver linear funcionou!")
                    except:
                        if self.rank == 0:
                            print(f"⚠️  Solver linear também falhou, mantendo solução anterior")
                elif self.rank == 0 and n_iterations > 10:
                    print(f"⚠️  Convergiu em {n_iterations} iterações (lento)")
                    
            except Exception as e:
                if self.rank == 0:
                    print(f"⚠️  Erro no sistema simples: {e}")
    
    def solve_exothermic_timestep(self, dt_val, current_time, bcs):
        """Resolve sistema exotérmico com duas equações acopladas"""
        dt = Constant(self.mesh, PETSc.ScalarType(dt_val))
        theta = Constant(self.mesh, PETSc.ScalarType(self.theta))
        
        # Medidas de integração
        dx_tags = dx(domain=self.mesh, subdomain_data=self.cell_tags)
        ds_tags = ds(domain=self.mesh, subdomain_data=self.facet_tags)
        
        try:
            # 1. RESOLVER TEMPERATURA (usando teq_n fixo)
            F_T = self.build_temperature_equation(dt, theta, dx_tags, ds_tags, current_time)
            problem_T = NonlinearProblem(F_T, self.T, bcs)
            solver_T = NewtonSolver(self.comm, problem_T)
            solver_T.convergence_criterion = "incremental"
            solver_T.rtol = 1e-6
            solver_T.max_it = 50
            
            n_iter_T, conv_T = solver_T.solve(self.T)
            
            # 2. RESOLVER TEMPO EQUIVALENTE (usando T novo)
            F_teq = self.build_teq_equation(dt, dx_tags)
            problem_teq = NonlinearProblem(F_teq, self.teq, [])
            solver_teq = NewtonSolver(self.comm, problem_teq)
            solver_teq.convergence_criterion = "incremental"
            solver_teq.rtol = 1e-6
            solver_teq.max_it = 20
            
            n_iter_teq, conv_teq = solver_teq.solve(self.teq)
            
            if self.rank == 0 and (not conv_T or not conv_teq):
                print(f"⚠️  Sistema exotérmico: T({n_iter_T})={'✓' if conv_T else '✗'}, teq({n_iter_teq})={'✓' if conv_teq else '✗'}")
                
        except Exception as e:
            if self.rank == 0:
                print(f"⚠️  Erro no sistema exotérmico: {e}")
    
    def build_temperature_equation(self, dt, theta, dx_tags, ds_tags, current_time):
        """Constrói equação da temperatura para sistema exotérmico"""
        F_T = 0
        
        # Processar apenas domínios ativos
        for domain_id in self.discovered_cell_tags:
            if domain_id not in self.active_layers or not self.active_layers[domain_id]:
                continue
            
            mat_props = self.get_material_properties(domain_id)
            
            rho = Constant(self.mesh, PETSc.ScalarType(mat_props['densidade']))
            k = Constant(self.mesh, PETSc.ScalarType(mat_props['condutividade']))
            cp = Constant(self.mesh, PETSc.ScalarType(mat_props['calor_especifico']))
            
            # Termo temporal
            F_T += rho * cp * (self.T - self.T_n) / dt * self.v_T * dx_tags(domain_id)
            
            # Termo de difusão
            T_theta = theta * self.T + (1 - theta) * self.T_n
            F_T += k * dot(grad(T_theta), grad(self.v_T)) * dx_tags(domain_id)
            
            # Termo exotérmico (se material gera calor)
            if mat_props['gera_calor']:
                # dQ/dteq * (teq - teq_n) / dt (usando valores anteriores)
                dQ_val = self.dQ_dteq_numeric(self.teq_n, mat_props['densidade'], mat_props['calor_especifico'])
                if dQ_val > 0:  # Evitar valores negativos
                    dQ = Constant(self.mesh, PETSc.ScalarType(dQ_val))
                    F_T -= dQ * (self.teq - self.teq_n) / dt * self.v_T * dx_tags(domain_id)
        
        # Aplicar APENAS condições de contorno dos contornos ATIVOS
        contornos_ativos = self.get_active_boundaries(current_time)
        
        for boundary_tag in contornos_ativos:
            if boundary_tag in self.boundary_conditions:
                bc_config = self.boundary_conditions[boundary_tag]
                if bc_config['tipo'] == 'conveccao':
                    h_val = bc_config['h']
                    T_ext_val = bc_config['t_ext']
                    
                    h = Constant(self.mesh, PETSc.ScalarType(h_val))
                    T_ext = Constant(self.mesh, PETSc.ScalarType(T_ext_val))
                    
                    T_boundary = theta * self.T + (1 - theta) * self.T_n
                    F_T += h * (T_boundary - T_ext) * self.v_T * ds_tags(boundary_tag)
        
        return F_T
    
    def build_teq_equation(self, dt, dx_tags):
        """Constrói equação do tempo equivalente"""
        F_teq = 0
        
        # Processar apenas domínios ativos
        for domain_id in self.discovered_cell_tags:
            if domain_id not in self.active_layers or not self.active_layers[domain_id]:
                continue
            
            mat_props = self.get_material_properties(domain_id)
            
            if mat_props['gera_calor']:
                # dteq/dt = exp(Ea/R * (1/T0 - 1/T))
                # (teq - teq_n) / dt = arrhenius(T_n)
                F_teq += (self.teq - self.teq_n) / dt * self.v_teq * dx_tags(domain_id)
                
                # Fator de Arrhenius (usando temperatura anterior para estabilidade)
                arrhenius_val = self.arrhenius_factor_numeric(self.T_n)
                arrhenius = Constant(self.mesh, PETSc.ScalarType(arrhenius_val))
                F_teq -= arrhenius * self.v_teq * dx_tags(domain_id)
            else:
                # Material sem geração: teq permanece constante
                F_teq += (self.teq - self.teq_n) * self.v_teq * dx_tags(domain_id)
        
        return F_teq
    
    def dQ_dteq_numeric(self, teq_func, rho, cp):
        """Versão numérica de dQ/dteq para usar em Constant"""
        # Usar valor médio do tempo equivalente anterior
        teq_avg = np.mean(teq_func.x.array[:]) if len(teq_func.x.array) > 0 else 0.0
        teq_avg = max(teq_avg, 1e-6)  # Evitar zero
        
        term1 = (self.Tad_inf1 * self.c1 * teq_avg**(self.c1-1) * self.a1**self.c1 / 
                (self.a1**self.c1 + teq_avg**self.c1)**2)
        term2 = (self.Tad_inf2 * self.c2 * teq_avg**(self.c2-1) * self.a2**self.c2 / 
                (self.a2**self.c2 + teq_avg**self.c2)**2)
        return rho * cp * (term1 + term2)
    
    def arrhenius_factor_numeric(self, T_func):
        """Versão numérica do fator de Arrhenius"""
        # Usar temperatura média
        T_avg = np.mean(T_func.x.array[:]) if len(T_func.x.array) > 0 else 20.0
        return np.exp(self.EaR * (1/298.15 - 1/(T_avg + 273.15)))
    
    def save_results(self, time_step, current_time):
        """Salva resultados"""
        try:
            if self.has_exothermic:
                # Salvar temperatura e tempo equivalente (separadamente)
                output_temp = f"{self.output_dir}/temperatura_{time_step:04d}.xdmf"
                output_teq = f"{self.output_dir}/tempo_equiv_{time_step:04d}.xdmf"
                
                # Salvar temperatura
                with io.XDMFFile(self.comm, output_temp, "w") as xdmf:
                    self.T.name = "Temperatura"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.T, current_time)
                
                # Salvar tempo equivalente
                with io.XDMFFile(self.comm, output_teq, "w") as xdmf:
                    self.teq.name = "TempoEquivalente"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.teq, current_time)
            else:
                # Sistema simples: apenas temperatura
                output_file = f"{self.output_dir}/temperatura_{time_step:04d}.xdmf"
                with io.XDMFFile(self.comm, output_file, "w") as xdmf:
                    self.T.name = "Temperatura"
                    xdmf.write_mesh(self.mesh)
                    xdmf.write_function(self.T, current_time)
        except Exception as e:
            if self.rank == 0:
                print(f"⚠️  Erro ao salvar: {e}")
    
    def run_simulation(self):
        """Executa simulação com sistema de duas fases (igual aos exemplos)"""
        if self.rank == 0:
            if self.has_exothermic:
                print("\n🚀 Iniciando simulação EXOTÉRMICA com duas fases...")
                print(f"🔥 Fase 1: dt={self.delta_t_refinado/3600:.1f}h até 48h")
                print(f"🔥 Fase 2: dt={self.delta_t/3600:.1f}h até {self.tempo_final/86400:.1f} dias")
            else:
                print("\n🚀 Iniciando simulação SIMPLES com etapas construtivas...")
        
        current_time = 0.0
        time_step = 0
        T_split = 48 * 3600  # 48 horas em segundos
        
        while current_time < self.tempo_final:
            # Sistema de duas fases (igual aos exemplos)
            if current_time < T_split:
                dt_val = self.delta_t_refinado  # dt pequeno até 48h
            else:
                dt_val = self.delta_t  # dt maior após 48h
            
            if current_time + dt_val > self.tempo_final:
                dt_val = self.tempo_final - current_time
            
            # Resolver
            self.solve_timestep(dt_val, current_time)
            
            # Atualizar soluções anteriores
            if self.has_exothermic:
                # Sistema exotérmico: atualizar Tp e teq
                self.T_n.x.array[:] = self.T.x.array[:]
                self.teq_n.x.array[:] = self.teq.x.array[:]
            else:
                # Sistema simples: atualizar apenas T
                self.Tn.x.array[:] = self.T.x.array[:]
            
            current_time += dt_val
            time_step += 1
            
            # Salvar resultados e log
            save_frequency = 1 if current_time < T_split else 5  # Mais frequente na fase 1
            
            if time_step % save_frequency == 0:
                self.save_results(time_step, current_time)
                
                if self.rank == 0:
                    if self.has_exothermic:
                        # Log para sistema exotérmico
                        print(f"⏰ Tempo: {current_time/3600:.1f}h, Step: {time_step}, dt: {dt_val/3600:.1f}h")
                    else:
                        # Log para sistema simples
                        progress = current_time / self.tempo_final * 100
                        print(f"📊 {progress:.1f}% - {current_time/3600:.1f}h")
        
        if self.rank == 0:
            print(f"\n✅ Simulação concluída!")
            print(f"⏱️  Tempo total: {current_time/86400:.1f} dias")
            print(f"📊 Total de passos: {time_step}")
            print(f"📁 Resultados salvos em: {self.output_dir}")


def main():
    """Função principal genérica"""
    yaml_file = "barragem1.yaml"
    xdmf_file = "malha_barragem.xdmf"
    
    if not os.path.exists(yaml_file) or not os.path.exists(xdmf_file):
        print("❌ Arquivos não encontrados")
        return
    
    solver = GenericThermalFEniCSx(yaml_file, xdmf_file)
    solver.run_simulation()


if __name__ == "__main__":
    main() 