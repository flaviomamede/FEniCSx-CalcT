#!/usr/bin/env python3
"""
Analisador Stagewise GENÉRICO para Barragem FEniCSx

NOVIDADES DESTA VERSÃO:
- Mapeamento dinâmico de Nomes para IDs de Physical Groups a partir do YAML.
- Suporte à desativação explícita de contornos via chave 'desativado_pela_camada'.
- Vetor de tempo refinado cirurgicamente em torno dos 'birth events'.
- DUPLA SAÍDA: JSON com numeração .msh (para conferência) e .xdmf/.h5 (para solver)

Este script extrai e analisa:
1. Vetor de tempo utilizado.
2. Physical Groups ativos por bloco de tempo.
3. Nós e elementos de domínio e contorno ativos por bloco de tempo.
"""

import numpy as np
import yaml
from mpi4py import MPI
from dolfinx import io
from pathlib import Path
import json
from datetime import datetime

class AnalisadorStagewiseGenerico:
    """
    Classe para análise completa do processo stage-wise com mapeamento dinâmico.
    """
    
    def __init__(self, yaml_file, xdmf_file, msh_file=None):
        self.yaml_file = yaml_file
        self.xdmf_file = xdmf_file
        self.msh_file = msh_file or xdmf_file.replace('.xdmf', '.msh')
        self.comm = MPI.COMM_WORLD
        
        # Carregar configurações do YAML
        with open(yaml_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        # Carregar a malha computacional
        self._carregar_malha()
        
        # Carregar correspondência .msh <-> .xdmf/.h5
        self._carregar_correspondencia_msh()
        
        ## NOVO: Dicionários que serão preenchidos dinamicamente
        self.physical_surfaces = {}
        self.physical_lines = {}
        
        ## NOVO: Chamar o método que constrói os mapeamentos
        self._mapear_nomes_para_ids()
        
        # Estruturas para armazenar os resultados da análise
        self.vetor_tempo = []
        self.blocos_tempo = []
        self.analise_resultados = {}
        
    def _carregar_malha(self):
        """
        Carrega a malha e os meshtags
        """
        print(f"Carregando malha: {self.xdmf_file}")
        with io.XDMFFile(self.comm, self.xdmf_file, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            self.mesh.topology.create_entities(1)
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        print(f"Malha carregada: {self.mesh.topology.index_map(2).size_local} células, "
              f"{self.mesh.topology.index_map(1).size_local} facetas")

    def _carregar_correspondencia_msh(self):
        """
        Carrega correspondência entre elementos .msh e .xdmf/.h5
        VERSÃO ROBUSTA: Suporta diferentes tipos de elementos e tolerâncias
        """
        print(f"Carregando correspondência do arquivo: {self.msh_file}")
        
        # Configurações robustas
        TOLERANCIA_COORDENADAS = 0.001  # Pode ser ajustado conforme necessário
        TIPOS_ELEMENTOS_SUPORTADOS = {
            2: 3,   # Triângulo - 3 nós
            3: 4,   # Quadrilátero - 4 nós
            4: 4,   # Tetraédrico - 4 nós
            5: 8,   # Hexaédrico - 8 nós
        }
        
        # Ler arquivo .msh
        try:
            with open(self.msh_file, 'r') as f:
                lines = f.readlines()
        except FileNotFoundError:
            print(f"❌ Arquivo .msh não encontrado: {self.msh_file}")
            self.correspondencia_msh_xdmf = {}
            self.correspondencia_xdmf_msh = {}
            return

        # Extrair elementos e nós do .msh
        in_elements = False
        in_nodes = False
        elementos_msh = []
        nodes_msh = {}

        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            
            if line == '$Nodes':
                in_nodes = True
                continue
            elif line == '$EndNodes':
                in_nodes = False
                continue
            elif line == '$Elements':
                in_elements = True
                continue
            elif line == '$EndElements':
                in_elements = False
                continue
            
            if in_nodes and line and not line.startswith('$'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        node_id = int(parts[0])
                        x, y = float(parts[1]), float(parts[2])
                        nodes_msh[node_id] = (x, y)
                    except (ValueError, IndexError) as e:
                        print(f"  Aviso: erro ao processar nó na linha {line_num}: {e}")
            
            if in_elements and line and not line.startswith('$'):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        elem_id = int(parts[0])
                        elem_type = int(parts[1])
                        physical_group = int(parts[3])
                        
                        # Suporte robusto para diferentes tipos de elementos
                        if elem_type in TIPOS_ELEMENTOS_SUPORTADOS:
                            num_nodes = TIPOS_ELEMENTOS_SUPORTADOS[elem_type]
                            if len(parts) >= 5 + num_nodes:
                                node_ids = [int(parts[j]) for j in range(5, 5 + num_nodes)]
                                elementos_msh.append((elem_id, elem_type, physical_group, node_ids))
                        else:
                            print(f"  Aviso: tipo de elemento {elem_type} não suportado (linha {line_num})")
                    except (ValueError, IndexError) as e:
                        print(f"  Aviso: erro ao processar elemento na linha {line_num}: {e}")

        print(f"  📊 Elementos .msh processados: {len(elementos_msh)}")
        print(f"  📊 Nós .msh processados: {len(nodes_msh)}")

        # Calcular centroides dos elementos do .msh
        elementos_msh_centroids = []
        for elem_id, elem_type, pg, node_ids in elementos_msh:
            if all(nid in nodes_msh for nid in node_ids):
                coords = [nodes_msh[nid] for nid in node_ids]
                centroid_x = sum(coord[0] for coord in coords) / len(coords)
                centroid_y = sum(coord[1] for coord in coords) / len(coords)
                elementos_msh_centroids.append((elem_id, elem_type, pg, centroid_x, centroid_y))
            else:
                print(f"  Aviso: nós não encontrados para elemento {elem_id}")

        # Obter centroides dos elementos do .xdmf/.h5
        coords = self.mesh.geometry.x
        cell_to_vertex = self.mesh.topology.connectivity(self.mesh.topology.dim, 0)
        
        elementos_xdmf_centroids = []
        for elem_id in range(self.mesh.topology.index_map(2).size_local):
            nodes = cell_to_vertex.links(elem_id)
            centroid = np.mean([coords[node] for node in nodes], axis=0)
            physical_id = self.cell_tags.values[elem_id]
            elementos_xdmf_centroids.append((elem_id, physical_id, centroid[0], centroid[1]))

        print(f"  📊 Elementos .xdmf/.h5 processados: {len(elementos_xdmf_centroids)}")

        # Criar correspondência baseada nos centroides
        self.correspondencia_msh_xdmf = {}
        self.correspondencia_xdmf_msh = {}
        correspondencias_problematicas = []
        
        for msh_id, msh_type, msh_pg, msh_x, msh_y in elementos_msh_centroids:
            correspondencia_encontrada = False
            
            for xdmf_id, xdmf_pg, xdmf_x, xdmf_y in elementos_xdmf_centroids:
                distancia = np.sqrt((msh_x - xdmf_x)**2 + (msh_y - xdmf_y)**2)
                
                if distancia < TOLERANCIA_COORDENADAS and msh_pg == xdmf_pg:
                    self.correspondencia_msh_xdmf[msh_id] = xdmf_id
                    self.correspondencia_xdmf_msh[xdmf_id] = msh_id
                    correspondencia_encontrada = True
                    break
            
            if not correspondencia_encontrada:
                correspondencias_problematicas.append((msh_id, msh_type, msh_pg, msh_x, msh_y))
        
        print(f"  ✅ Correspondências criadas: {len(self.correspondencia_msh_xdmf)}")
        
        if correspondencias_problematicas:
            print(f"  ⚠️  Elementos sem correspondência: {len(correspondencias_problematicas)}")
            for elem_id, elem_type, pg, x, y in correspondencias_problematicas[:3]:
                print(f"    Elemento {elem_id} (tipo {elem_type}, PG {pg}) em ({x:.6f}, {y:.6f})")
        
        # Verificar se a correspondência foi bem-sucedida
        if len(self.correspondencia_msh_xdmf) == 0:
            print("  ❌ ERRO: Nenhuma correspondência encontrada!")
            print("  💡 Possíveis causas:")
            print("     - Arquivos .msh e .xdmf/.h5 não correspondem")
            print("     - Tolerância muito baixa")
            print("     - Physical Groups diferentes")
        else:
            taxa_sucesso = len(self.correspondencia_msh_xdmf) / len(elementos_msh_centroids) * 100
            print(f"  📊 Taxa de sucesso: {taxa_sucesso:.1f}%")
            
            if taxa_sucesso < 95:
                print("  ⚠️  Taxa de sucesso baixa - verifique os arquivos!")

    def _mapear_nomes_para_ids(self):
        """
        Constrói os dicionários de mapeamento (nome -> ID) lendo do próprio YAML.
        """
        print("\nConstruindo mapeamentos dinâmicos de Nomes para IDs...")
        
        # Mapeamento para Physical Surfaces (Domínios)
        for camada_mat in self.config.get('camadas_material', []):
            nome = camada_mat.get('nome')
            # Extrai o ID numérico do final do nome (ex: "camada_material_10" -> 10)
            try:
                id_numerico = int(nome.split('_')[-1])
                self.physical_surfaces[nome] = id_numerico
            except (ValueError, IndexError):
                print(f"  Aviso: não foi possível extrair ID do nome '{nome}'.")
        
        # Mapeamento para Physical Lines (Contornos)
        for contorno in self.config.get('contornos', []):
            nome = contorno.get('nome')
            contorno_id = contorno.get('id')
            if nome and contorno_id is not None:
                self.physical_lines[nome] = contorno_id
        
        print(f"  📋 {len(self.physical_surfaces)} mapeamentos de superfície criados.")
        print(f"  📋 {len(self.physical_lines)} mapeamentos de linha criados.")
        
        # Mostrar mapeamentos criados para debug
        if self.physical_surfaces:
            print("  🔹 Superfícies mapeadas:", dict(list(self.physical_surfaces.items())[:5]))
        if self.physical_lines:
            print("  🔹 Linhas mapeadas:", dict(list(self.physical_lines.items())[:5]))

    def gerar_vetor_tempo(self):
        """
        Gera o vetor de tempo baseado na configuração YAML, com passos refinados
        em torno dos eventos de 'birth' de cada camada.
        """
        print("\n=== 1. GERANDO VETOR DE TEMPO (LÓGICA REFINADA) ===")
        
        general = self.config['general']
        tempo_final = general['tempo_final']
        delta_t = general['delta_t']
        delta_t_refinado = general['delta_t_refinado']
        
        # 1. Coletar todos os pontos de tempo "padrão"
        pontos_padrao = np.arange(0, tempo_final + delta_t, delta_t)

        # 2. Coletar todos os pontos "especiais" (críticos e refinados)
        pontos_especiais = {0.0, tempo_final}

        for camada in self.config['camadas']:
            birth_time = camada['birth']
            pontos_especiais.add(birth_time)
            if birth_time > 0:
                pontos_especiais.add(birth_time - delta_t_refinado)
            pontos_especiais.add(birth_time + delta_t_refinado)

        # 3. Combinar tudo, remover duplicatas e ordenar
        todos_os_pontos = set(pontos_padrao)
        todos_os_pontos.update(pontos_especiais)
        
        self.vetor_tempo = np.array(sorted([p for p in todos_os_pontos if p <= tempo_final]))

        print(f"Vetor de tempo gerado com {len(self.vetor_tempo)} pontos.")
        
        # Mostrar refinamentos em torno dos births
        birth_times = [camada['birth'] for camada in self.config['camadas']]
        print(f"Tempos de birth das camadas: {birth_times}")
        
        return self.vetor_tempo
    
    def identificar_blocos_tempo(self):
        """
        Identifica os blocos de tempo baseado no birth das camadas
        """
        print("\n=== 2. IDENTIFICANDO BLOCOS DE TEMPO ===")
        
        birth_times = sorted(set([c['birth'] for c in self.config['camadas']]))
        
        tempo_final = self.config['general']['tempo_final']
        if tempo_final not in birth_times:
            birth_times.append(tempo_final)
        
        self.blocos_tempo = []
        for i in range(len(birth_times) - 1):
            bloco = {
                'id': i + 1,
                'inicio': birth_times[i],
                'fim': birth_times[i + 1],
                'duracao': birth_times[i + 1] - birth_times[i]
            }
            self.blocos_tempo.append(bloco)
        
        print(f"Identificados {len(self.blocos_tempo)} blocos de tempo.")
        for bloco in self.blocos_tempo:
            print(f"  Bloco {bloco['id']}: {bloco['inicio']} - {bloco['fim']} s "
                  f"({bloco['duracao']/86400:.1f} dias)")
        
        return self.blocos_tempo
    
    def obter_camadas_ativas(self, tempo):
        """
        Obtém as camadas ativas em um dado tempo
        """
        camadas_ativas = []
        for camada in self.config['camadas']:
            if camada['birth'] <= tempo and (camada['death'] is None or tempo < camada['death']):
                camadas_ativas.append(camada['nome'])
        return camadas_ativas
    
    def obter_physical_groups_ativos(self, camadas_ativas):
        """
        Obtém os physical groups ativos baseado nas camadas ativas
        COM SUPORTE À DESATIVAÇÃO EXPLÍCITA e MAPEAMENTO DINÂMICO
        """
        surfaces_ativas = []
        lines_ativas = []
        
        # Domínios (surfaces) usando mapeamento dinâmico
        for camada_material in self.config['camadas_material']:
            if camada_material['camada'] in camadas_ativas:
                nome = camada_material['nome']
                if nome in self.physical_surfaces:
                    surfaces_ativas.append(self.physical_surfaces[nome])
        
        # Contornos (lines) usando mapeamento dinâmico
        for contorno in self.config['contornos']:
            nasce_com = contorno.get('nasce_com_camada')
            desativado_por = contorno.get('desativado_pela_camada')
            
            # Verificar se o contorno deve estar ativo
            if nasce_com in camadas_ativas and (not desativado_por or desativado_por not in camadas_ativas):
                nome = contorno.get('nome')
                if nome in self.physical_lines:
                    lines_ativas.append(self.physical_lines[nome])
        
        return {'surfaces': sorted(surfaces_ativas), 'lines': sorted(lines_ativas)}
    
    def obter_elementos_e_nos(self, surface_ids=None, line_ids=None):
        """
        Obtém elementos e nós associados aos physical groups
        """
        resultado = {
            'elementos_dominio': [],
            'nos_dominio': set(),
            'elementos_contorno': [],
            'nos_contorno': set()
        }
        
        if surface_ids:
            for sid in surface_ids:
                elementos = self.cell_tags.find(sid)
                resultado['elementos_dominio'].extend(elementos)
                for elem in elementos:
                    resultado['nos_dominio'].update(self.mesh.topology.connectivity(2, 0).links(elem))
        
        if line_ids:
            for lid in line_ids:
                elementos = self.facet_tags.find(lid)
                resultado['elementos_contorno'].extend(elementos)
                for elem in elementos:
                    resultado['nos_contorno'].update(self.mesh.topology.connectivity(1, 0).links(elem))
        
        # Converter sets para listas ordenadas
        for key in ['elementos_dominio', 'nos_dominio', 'elementos_contorno', 'nos_contorno']:
            if isinstance(resultado[key], set):
                resultado[key] = sorted(list(resultado[key]))
            else:
                resultado[key] = sorted(resultado[key])
        
        return resultado
    
    def calcular_diferencas(self, atual, anterior):
        """
        Calcula as diferenças entre dois conjuntos (entradas e saídas)
        """
        atual_set, anterior_set = set(atual), set(anterior)
        return {
            'entradas': sorted(list(atual_set - anterior_set)),
            'saidas': sorted(list(anterior_set - atual_set))
        }
    
    def analisar_todos_blocos(self):
        """
        Análise completa de todos os blocos de tempo
        """
        print("\n=== 3. ANALISANDO TODOS OS BLOCOS DE TEMPO ===")
        
        for i, bloco in enumerate(self.blocos_tempo):
            print(f"\nAnalisando Bloco {bloco['id']}...")
            
            tempo_representativo = bloco['inicio']
            camadas_ativas = self.obter_camadas_ativas(tempo_representativo)
            physical_groups = self.obter_physical_groups_ativos(camadas_ativas)
            elementos_nos = self.obter_elementos_e_nos(
                surface_ids=physical_groups['surfaces'],
                line_ids=physical_groups['lines']
            )
            
            # Calcular diferenças
            diferencas = {}
            if i > 0:
                bloco_anterior = self.analise_resultados[f'bloco_{i}']
                
                diferencas['surfaces'] = self.calcular_diferencas(
                    physical_groups['surfaces'],
                    bloco_anterior['physical_groups']['surfaces']
                )
                
                diferencas['lines'] = self.calcular_diferencas(
                    physical_groups['lines'],
                    bloco_anterior['physical_groups']['lines']
                )
                
                diferencas['elementos_dominio'] = self.calcular_diferencas(
                    elementos_nos['elementos_dominio'],
                    bloco_anterior['elementos_nos']['elementos_dominio']
                )
                
                diferencas['nos_dominio'] = self.calcular_diferencas(
                    elementos_nos['nos_dominio'],
                    bloco_anterior['elementos_nos']['nos_dominio']
                )
                
                diferencas['elementos_contorno'] = self.calcular_diferencas(
                    elementos_nos['elementos_contorno'],
                    bloco_anterior['elementos_nos']['elementos_contorno']
                )
                
                diferencas['nos_contorno'] = self.calcular_diferencas(
                    elementos_nos['nos_contorno'],
                    bloco_anterior['elementos_nos']['nos_contorno']
                )
            
            # Armazenar resultados
            self.analise_resultados[f'bloco_{i+1}'] = {
                'info_bloco': bloco,
                'camadas_ativas': camadas_ativas,
                'physical_groups': physical_groups,
                'elementos_nos': elementos_nos,
                'diferencas': diferencas
            }
    
    def gerar_relatorio(self):
        """
        Gera relatório detalhado da análise
        """
        print("\n" + "="*80)
        print("RELATÓRIO COMPLETO DA ANÁLISE STAGE-WISE GENÉRICA")
        print("="*80)
        
        # 1. Informações gerais
        print(f"\n1. INFORMAÇÕES GERAIS")
        print(f"   Arquivo YAML: {self.yaml_file}")
        print(f"   Arquivo XDMF: {self.xdmf_file}")
        print(f"   Mapeamentos criados: {len(self.physical_surfaces)} superfícies, {len(self.physical_lines)} linhas")
        print(f"   Vetor de tempo: {len(self.vetor_tempo)} pontos")
        
        # 2. Análise por bloco
        for i, (bloco_key, bloco_data) in enumerate(self.analise_resultados.items()):
            bloco_info = bloco_data['info_bloco']
            print(f"\n2.{i+1}. BLOCO {bloco_info['id']} ({bloco_info['inicio']} - {bloco_info['fim']} s)")
            print(f"     Duração: {bloco_info['duracao']/86400:.1f} dias")
            print(f"     Camadas ativas: {bloco_data['camadas_ativas']}")
            
            pg = bloco_data['physical_groups']
            print(f"     Physical Surfaces (domínios): {pg['surfaces']}")
            print(f"     Physical Lines (contornos): {pg['lines']}")
            
            en = bloco_data['elementos_nos']
            print(f"     Elementos: {len(en['elementos_dominio'])} domínio, {len(en['elementos_contorno'])} contorno")
            print(f"     Nós: {len(en['nos_dominio'])} domínio, {len(en['nos_contorno'])} contorno")
            
            # Mostrar mudanças
            if bloco_data['diferencas']:
                diff = bloco_data['diferencas']
                if diff['surfaces']['entradas'] or diff['surfaces']['saidas']:
                    print(f"     🔄 Surfaces: +{diff['surfaces']['entradas']} -{diff['surfaces']['saidas']}")
                if diff['lines']['entradas'] or diff['lines']['saidas']:
                    print(f"     🔄 Lines: +{diff['lines']['entradas']} -{diff['lines']['saidas']}")
        
        print(f"\n{'='*80}")
        print(f"Análise GENÉRICA concluída em {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*80}")
    
    def converter_elementos_xdmf_para_msh(self, elementos_xdmf):
        """
        Converte lista de elementos da numeração XDMF para MSH
        """
        elementos_msh = []
        for elem_xdmf in elementos_xdmf:
            if elem_xdmf in self.correspondencia_xdmf_msh:
                elementos_msh.append(self.correspondencia_xdmf_msh[elem_xdmf])
        return sorted(elementos_msh)

    def criar_resultados_versao_msh(self):
        """
        Cria versão dos resultados com numeração MSH
        """
        analise_resultados_msh = {}
        
        for bloco_key, bloco_data in self.analise_resultados.items():
            # Converter elementos de domínio
            elementos_dominio_xdmf = bloco_data['elementos_nos']['elementos_dominio']
            elementos_dominio_msh = self.converter_elementos_xdmf_para_msh(elementos_dominio_xdmf)
            
            # Converter elementos de contorno
            elementos_contorno_xdmf = bloco_data['elementos_nos']['elementos_contorno']
            elementos_contorno_msh = self.converter_elementos_xdmf_para_msh(elementos_contorno_xdmf)
            
            # Criar nova estrutura com numeração MSH
            analise_resultados_msh[bloco_key] = {
                'info_bloco': bloco_data['info_bloco'],
                'camadas_ativas': bloco_data['camadas_ativas'],
                'physical_groups': bloco_data['physical_groups'],
                'elementos_nos': {
                    'elementos_dominio': elementos_dominio_msh,
                    'nos_dominio': bloco_data['elementos_nos']['nos_dominio'],  # Nós mantém numeração
                    'elementos_contorno': elementos_contorno_msh,
                    'nos_contorno': bloco_data['elementos_nos']['nos_contorno']  # Nós mantém numeração
                },
                'diferencas': bloco_data['diferencas']
            }
        
        return analise_resultados_msh

    def salvar_resultados_duplo(self):
        """
        Salva os resultados em dois arquivos JSON: um com numeração MSH e outro com XDMF
        """
        def converter_para_json(obj):
            if isinstance(obj, np.generic):
                return obj.item()
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        # === ARQUIVO 1: NUMERAÇÃO XDMF (PARA O SOLVER) ===
        resultados_xdmf = {
            'info_geral': {
                'yaml_file': self.yaml_file,
                'xdmf_file': self.xdmf_file,
                'msh_file': self.msh_file,
                'tipo_numeracao': 'XDMF/H5 (para solver)',
                'mapeamentos': {
                    'physical_surfaces': self.physical_surfaces,
                    'physical_lines': self.physical_lines
                },
                'correspondencia_disponivel': len(self.correspondencia_msh_xdmf)
            },
            'vetor_tempo': self.vetor_tempo.tolist(),
            'blocos_tempo': self.blocos_tempo,
            'analise_resultados': self.analise_resultados
        }
        
        arquivo_xdmf = "analise_stagewise_xdmf.json"
        with open(arquivo_xdmf, 'w', encoding='utf-8') as f:
            json.dump(resultados_xdmf, f, indent=2, ensure_ascii=False, default=converter_para_json)
        
        # === ARQUIVO 2: NUMERAÇÃO MSH (PARA CONFERÊNCIA) ===
        analise_resultados_msh = self.criar_resultados_versao_msh()
        
        resultados_msh = {
            'info_geral': {
                'yaml_file': self.yaml_file,
                'xdmf_file': self.xdmf_file,
                'msh_file': self.msh_file,
                'tipo_numeracao': 'MSH (para conferência)',
                'mapeamentos': {
                    'physical_surfaces': self.physical_surfaces,
                    'physical_lines': self.physical_lines
                },
                'correspondencia_disponivel': len(self.correspondencia_msh_xdmf)
            },
            'vetor_tempo': self.vetor_tempo.tolist(),
            'blocos_tempo': self.blocos_tempo,
            'analise_resultados': analise_resultados_msh,
            'correspondencia_msh_xdmf': self.correspondencia_msh_xdmf
        }
        
        arquivo_msh = "analise_stagewise_msh.json"
        with open(arquivo_msh, 'w', encoding='utf-8') as f:
            json.dump(resultados_msh, f, indent=2, ensure_ascii=False, default=converter_para_json)
        
        print(f"\n📄 Resultados salvos em:")
        print(f"  🔧 Para o SOLVER: {arquivo_xdmf}")
        print(f"  ✅ Para CONFERÊNCIA: {arquivo_msh}")
        
        return arquivo_xdmf, arquivo_msh
    
    def executar_analise_completa(self):
        """
        Executa a análise completa
        """
        print("INICIANDO ANÁLISE COMPLETA STAGE-WISE GENÉRICA")
        print("="*55)
        
        try:
            self.gerar_vetor_tempo()
            self.identificar_blocos_tempo()
            self.analisar_todos_blocos()
            self.gerar_relatorio()
            arquivo_xdmf, arquivo_msh = self.salvar_resultados_duplo()
            
            print("\n✅ ANÁLISE GENÉRICA CONCLUÍDA COM SUCESSO!")
            print(f"\n📋 RESUMO DOS ARQUIVOS GERADOS:")
            print(f"  🔧 {arquivo_xdmf} - Numeração XDMF/H5 (para usar no solver)")
            print(f"  ✅ {arquivo_msh} - Numeração MSH (para você conferir)")
            
        except Exception as e:
            print(f"\n❌ ERRO DURANTE A ANÁLISE: {e}")
            raise


def main():
    """
    Função principal
    """
    import sys
    
    # Verificar argumentos de linha de comando
    if len(sys.argv) != 3:
        print("❌ USO: python analisador_stagewise_generico.py <arquivo.xdmf> <arquivo.yaml>")
        print("    Exemplo: python analisador_stagewise_generico.py barragem2.xdmf barragem2_melhorado.yaml")
        return
    
    # Arquivos de entrada
    xdmf_file = sys.argv[1]
    yaml_file = sys.argv[2]
    
    # Verificar se os arquivos existem
    if not Path(xdmf_file).exists():
        print(f"❌ Arquivo XDMF não encontrado: {xdmf_file}")
        return
    
    if not Path(yaml_file).exists():
        print(f"❌ Arquivo YAML não encontrado: {yaml_file}")
        return
    
    # Criar analisador e executar
    analisador = AnalisadorStagewiseGenerico(yaml_file, xdmf_file)
    analisador.executar_analise_completa()


if __name__ == "__main__":
    main() 