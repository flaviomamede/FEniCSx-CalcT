"""
Exemplo complementar: Configuração de malha com tags para construção em camadas
Integração com GMSH para definição de camadas e condições de contorno
"""

import numpy as np
import dolfinx
import dolfinx.io
import dolfinx.mesh
from mpi4py import MPI
import gmsh

def criar_malha_barragem_camadas():
    """
    Cria malha de barragem com tags para construção em camadas usando GMSH
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    
    # Criar modelo
    gmsh.model.add("barragem_camadas")
    
    # Parâmetros da geometria
    largura_base = 100.0
    altura_total = 60.0
    num_camadas = 3
    altura_camada = altura_total / num_camadas
    
    # Pontos da geometria
    pontos = []
    
    # Base da barragem
    p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, 2.0)
    p2 = gmsh.model.geo.addPoint(largura_base, 0.0, 0.0, 2.0)
    pontos.extend([p1, p2])
    
    # Pontos das camadas
    for i in range(num_camadas + 1):
        y = i * altura_camada
        largura_nivel = largura_base - 2 * i * 5.0  # Formato trapezoidal
        x_esq = (largura_base - largura_nivel) / 2
        x_dir = x_esq + largura_nivel
        
        if i > 0:  # Não repetir pontos da base
            p_esq = gmsh.model.geo.addPoint(x_esq, y, 0.0, 2.0)
            p_dir = gmsh.model.geo.addPoint(x_dir, y, 0.0, 2.0)
            pontos.extend([p_esq, p_dir])
    
    # Criar linhas e definir camadas
    camadas_surfaces = []
    
    for i in range(num_camadas):
        # Índices dos pontos para esta camada
        if i == 0:
            # Primeira camada (base)
            p_base_esq = p1
            p_base_dir = p2
            p_topo_esq = pontos[2 + i * 2]
            p_topo_dir = pontos[3 + i * 2]
        else:
            # Camadas superiores
            p_base_esq = pontos[2 + (i-1) * 2]
            p_base_dir = pontos[3 + (i-1) * 2]
            p_topo_esq = pontos[2 + i * 2]
            p_topo_dir = pontos[3 + i * 2]
        
        # Criar linhas da camada
        l1 = gmsh.model.geo.addLine(p_base_esq, p_base_dir)
        l2 = gmsh.model.geo.addLine(p_base_dir, p_topo_dir)
        l3 = gmsh.model.geo.addLine(p_topo_dir, p_topo_esq)
        l4 = gmsh.model.geo.addLine(p_topo_esq, p_base_esq)
        
        # Criar loop da camada
        loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        
        # Criar superfície da camada
        surface = gmsh.model.geo.addPlaneSurface([loop])
        camadas_surfaces.append(surface)
        
        # Adicionar tag física para a camada
        gmsh.model.addPhysicalGroup(2, [surface], i + 1)
        gmsh.model.setPhysicalName(2, i + 1, f"camada_{i}")
    
    # Definir tags físicas para contornos
    # Tag 1: Base (fixa)
    gmsh.model.addPhysicalGroup(1, [1], 1)  # Linha da base
    gmsh.model.setPhysicalName(1, 1, "base_fixa")
    
    # Tags para contornos laterais de cada camada
    tag_contorno = 2
    for i in range(num_camadas):
        # Contornos laterais e superior de cada camada
        if i == 0:
            # Primeira camada: laterais
            gmsh.model.addPhysicalGroup(1, [2, 4], tag_contorno)
        else:
            # Camadas superiores: laterais
            linha_esq = 4 + (i-1) * 4 + 4
            linha_dir = 4 + (i-1) * 4 + 2
            gmsh.model.addPhysicalGroup(1, [linha_esq, linha_dir], tag_contorno)
        
        gmsh.model.setPhysicalName(1, tag_contorno, f"lateral_camada_{i}")
        tag_contorno += 1
    
    # Tag para topo (última camada)
    gmsh.model.addPhysicalGroup(1, [3 + (num_camadas-1) * 4], tag_contorno)
    gmsh.model.setPhysicalName(1, tag_contorno, "topo_final")
    
    # Sincronizar geometria
    gmsh.model.geo.synchronize()
    
    # Gerar malha
    gmsh.model.mesh.generate(2)
    
    # Salvar malha
    gmsh.write("barragem_camadas.msh")
    
    # Converter para dolfinx
    mesh, cell_tags, facet_tags = dolfinx.io.gmshio.read_from_msh(
        "barragem_camadas.msh", MPI.COMM_WORLD, gdim=2
    )
    
    gmsh.finalize()
    
    return mesh, cell_tags, facet_tags


def configurar_simulacao_barragem():
    """
    Configura simulação completa de barragem
    """
    # Criar malha
    mesh, cell_tags, facet_tags = criar_malha_barragem_camadas()
    
    # Parâmetros termicos do concreto
    parametros_termicos = {
        'condutividade': 2.5,      # W/m·K
        'densidade': 2400.0,       # kg/m³
        'calor_especifico': 1000.0, # J/kg·K
        'coef_expansao': 1e-5      # /K
    }
    
    # Configuração das camadas
    camadas_info = [
        {
            'id': 0,
            'tempo_lancamento': 0.0,
            'material': {
                'Q0': 800.0,      # Calor inicial [W/m³]
                'tau': 48.0,      # Constante de tempo [h]
                'E': 25e9,        # Módulo de elasticidade [Pa]
                'nu': 0.2,        # Coeficiente de Poisson
                'alpha': 1e-5     # Coeficiente de expansão térmica [/K]
            },
            'dirichlet': [{'tag': 1, 'valor': 20.0}],  # Base fixa
            'robin': [{'tag': 2, 'h': 12.0, 'T_amb': 15.0}]  # Lateral
        },
        {
            'id': 1,
            'tempo_lancamento': 168.0,  # 7 dias
            'material': {
                'Q0': 800.0,
                'tau': 48.0,
                'E': 25e9,
                'nu': 0.2,
                'alpha': 1e-5
            },
            'robin': [{'tag': 3, 'h': 12.0, 'T_amb': 15.0}]  # Lateral
        },
        {
            'id': 2,
            'tempo_lancamento': 336.0,  # 14 dias
            'material': {
                'Q0': 800.0,
                'tau': 48.0,
                'E': 25e9,
                'nu': 0.2,
                'alpha': 1e-5
            },
            'robin': [
                {'tag': 4, 'h': 12.0, 'T_amb': 15.0},  # Lateral
                {'tag': 5, 'h': 15.0, 'T_amb': 15.0}   # Topo
            ]
        }
    ]
    
    # Cronograma de construção
    cronograma = [
        {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},
        {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},
        {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2},
        {'tempo': 720.0, 'tipo': 'fim_simulacao'}  # 30 dias
    ]
    
    return mesh, cell_tags, facet_tags, parametros_termicos, camadas_info, cronograma


class ConstrucaoCamadasMelhorada:
    """
    Versão melhorada da classe com integração GMSH
    """
    
    def __init__(self, mesh, cell_tags, facet_tags, camadas_info, parametros_termicos):
        self.mesh = mesh
        self.cell_tags = cell_tags
        self.facet_tags = facet_tags
        self.camadas_info = camadas_info
        self.parametros_termicos = parametros_termicos
        
        # Configurar espaço funcional
        self.V = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))
        
        # Campos de solução
        self.T = dolfinx.fem.Function(self.V)
        self.T_old = dolfinx.fem.Function(self.V)
        
        # Campos de ativação por camada
        self.campos_ativacao = {}
        for camada in camadas_info:
            self.campos_ativacao[camada['id']] = dolfinx.fem.Function(self.V)
        
        # Campo de geração total
        self.campo_geracao = dolfinx.fem.Function(self.V)
        
        # Medidas de integração
        self.dx = ufl.Measure("dx", domain=mesh, subdomain_data=cell_tags)
        self.ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)
        
        # Estado da simulação
        self.tempo_atual = 0.0
        self.camada_atual = -1
        
        # Inicializar temperatura
        self.T.x.array[:] = 20.0  # Temperatura inicial
        self.T_old.x.array[:] = 20.0
    
    def obter_cells_camada(self, num_camada):
        """
        Retorna células da camada usando cell_tags
        """
        # A tag da camada é num_camada + 1 (pois tags começam em 1)
        tag_camada = num_camada + 1
        
        cells = []
        for i, tag in enumerate(self.cell_tags.values):
            if tag == tag_camada:
                cells.append(i)
        
        return np.array(cells)
    
    def ativar_camada(self, num_camada):
        """
        Ativa uma camada específica
        """
        print(f"Ativando camada {num_camada}")
        
        # Obter células da camada
        cells_camada = self.obter_cells_camada(num_camada)
        
        if len(cells_camada) == 0:
            print(f"Aviso: Nenhuma célula encontrada para camada {num_camada}")
            return
        
        # Ativar campo de ativação
        campo_ativacao = self.campos_ativacao[num_camada]
        campo_ativacao.x.array[:] = 0.0  # Resetar
        
        # Definir ativação apenas nas células da camada
        for cell in cells_camada:
            dofs_cell = self.V.dofmap.cell_dofs(cell)
            campo_ativacao.x.array[dofs_cell] = 1.0
        
        campo_ativacao.x.scatter_forward()
        
        # Atualizar campo de geração
        self.atualizar_geracao_calor()
        
        self.camada_atual = max(self.camada_atual, num_camada)
    
    def atualizar_geracao_calor(self):
        """
        Atualiza campo de geração de calor para todas as camadas ativas
        """
        self.campo_geracao.x.array[:] = 0.0
        
        for camada_info in self.camadas_info:
            camada_id = camada_info['id']
            
            if camada_id > self.camada_atual:
                continue
            
            # Calcular idade da camada
            idade = self.tempo_atual - camada_info['tempo_lancamento']
            
            if idade <= 0:
                continue
            
            # Calcular geração de calor
            Q0 = camada_info['material']['Q0']
            tau = camada_info['material']['tau']
            geracao = Q0 * np.exp(-idade / tau)
            
            # Aplicar geração nas células ativas da camada
            campo_ativacao = self.campos_ativacao[camada_id]
            geracao_camada = campo_ativacao.x.array * geracao
            
            # Somar ao campo total
            self.campo_geracao.x.array[:] += geracao_camada
        
        self.campo_geracao.x.scatter_forward()
    
    def definir_propriedades_efetivas(self):
        """
        Define propriedades efetivas baseadas nas camadas ativas
        """
        # Campos de propriedades
        k_eff = dolfinx.fem.Function(self.V)
        rho_cp_eff = dolfinx.fem.Function(self.V)
        
        k_eff.x.array[:] = 0.0
        rho_cp_eff.x.array[:] = 0.0
        
        # Somar contribuições de todas as camadas ativas
        for camada_info in self.camadas_info:
            camada_id = camada_info['id']
            
            if camada_id > self.camada_atual:
                continue
            
            # Propriedades da camada
            k = self.parametros_termicos['condutividade']
            rho = self.parametros_termicos['densidade']
            cp = self.parametros_termicos['calor_especifico']
            
            # Aplicar nas células ativas
            campo_ativacao = self.campos_ativacao[camada_id]
            
            k_eff.x.array[:] += campo_ativacao.x.array * k
            rho_cp_eff.x.array[:] += campo_ativacao.x.array * rho * cp
        
        k_eff.x.scatter_forward()
        rho_cp_eff.x.scatter_forward()
        
        return k_eff, rho_cp_eff
    
    def salvar_debug(self, nome_arquivo):
        """
        Salva campos para debug
        """
        arquivo = dolfinx.io.XDMFFile(MPI.COMM_WORLD, nome_arquivo, "w")
        arquivo.write_mesh(self.mesh)
        
        # Salvar temperatura
        arquivo.write_function(self.T, self.tempo_atual)
        
        # Salvar campo de geração
        arquivo.write_function(self.campo_geracao, self.tempo_atual)
        
        # Salvar campos de ativação
        for camada_id, campo in self.campos_ativacao.items():
            campo.name = f"ativacao_camada_{camada_id}"
            arquivo.write_function(campo, self.tempo_atual)
        
        arquivo.close()


if __name__ == "__main__":
    print("Configurando simulação de barragem com camadas...")
    
    # Configurar simulação
    mesh, cell_tags, facet_tags, parametros, camadas, cronograma = configurar_simulacao_barragem()
    
    # Criar simulação
    simulacao = ConstrucaoCamadasMelhorada(
        mesh, cell_tags, facet_tags, camadas, parametros
    )
    
    # Teste de ativação
    simulacao.ativar_camada(0)
    simulacao.tempo_atual = 100.0
    simulacao.atualizar_geracao_calor()
    
    # Salvar debug
    simulacao.salvar_debug("debug_barragem.xdmf")
    
    print("Configuração concluída!")