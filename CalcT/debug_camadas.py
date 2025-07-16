#!/usr/bin/env python3
"""
Analisador Stagewise GENÉRICO para Barragem FEniCSx

VERSÃO FINAL CORRIGIDA (2024-07-15)
- CORREÇÃO na leitura de MeshTags para mapeamento robusto de elementos 2D e 1D.
- Mapeamento e conversão de IDs para NÓS e ELEMENTOS.
- Mapeamento dinâmico de Nomes para IDs de Physical Groups a partir do YAML.
- DUPLA SAÍDA: JSON com numeração .msh (para conferência) e .xdmf/.h5 (para solver).
"""

import numpy as np
import yaml
from mpi4py import MPI
from dolfinx import io
from pathlib import Path
import json
from datetime import datetime
import sys

class AnalisadorStagewiseGenerico:
    """
    Classe para análise completa do processo stage-wise com mapeamento dinâmico.
    """
    
    def __init__(self, yaml_file, xdmf_file, msh_file=None):
        self.yaml_file = yaml_file
        self.xdmf_file = xdmf_file
        self.msh_file = msh_file or xdmf_file.replace('.xdmf', '.msh')
        self.comm = MPI.COMM_WORLD
        
        with open(yaml_file, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        self._carregar_malha()
        
        self.node_map_msh_to_xdmf, self.node_map_xdmf_to_msh = {}, {}
        self.elem2d_map_msh_to_xdmf, self.elem2d_map_xdmf_to_msh = {}, {}
        self.elem1d_map_msh_to_xdmf, self.elem1d_map_xdmf_to_msh = {}, {}
        
        self._carregar_correspondencia_msh()
        
        self.physical_surfaces, self.physical_lines = {}, {}
        self._mapear_nomes_para_ids()
        
        self.vetor_tempo, self.blocos_tempo, self.analise_resultados = [], [], {}

    def _carregar_malha(self):
        print(f"Carregando malha: {self.xdmf_file}")
        with io.XDMFFile(self.comm, self.xdmf_file, "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            self.mesh.topology.create_entities(1)
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        print(f"Malha carregada: {self.mesh.topology.index_map(self.mesh.topology.dim).size_global} células, "
              f"{self.mesh.topology.index_map(1).size_global} facetas")

    def _carregar_correspondencia_msh(self):
        print(f"Carregando correspondência do arquivo: {self.msh_file}")
        TOLERANCIA_COORDENADAS = 1e-6
        TIPOS_ELEMENTOS_SUPORTADOS = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 8: 3, 9: 6}
        
        with open(self.msh_file, 'r') as f: lines = f.readlines()

        nodes_msh, elementos2d_msh, elementos1d_msh = {}, [], []
        in_section = ''
        for line in lines:
            if '$Nodes' in line: in_section = 'nodes'; continue
            if '$Elements' in line: in_section = 'elements'; continue
            if '$End' in line: in_section = ''; continue
            if not in_section: continue
            
            parts = line.split()
            try:
                if in_section == 'nodes' and len(parts) >= 4:
                    nodes_msh[int(parts[0])] = tuple(float(p) for p in parts[1:4])
                elif in_section == 'elements' and len(parts) > 5:
                    elem_id, elem_type, num_tags = int(parts[0]), int(parts[1]), int(parts[2])
                    if elem_type in TIPOS_ELEMENTOS_SUPORTADOS:
                        pg = int(parts[3])
                        node_ids = [int(p) for p in parts[3 + num_tags:]]
                        if elem_type in [1, 8]:  # 1D
                            elementos1d_msh.append((elem_id, pg, node_ids))
                        else:  # 2D/3D
                            elementos2d_msh.append((elem_id, pg, node_ids))
            except (ValueError, IndexError): continue

        print(f"  📊 MSH Lido: {len(nodes_msh)} Nós, {len(elementos2d_msh)} Elem(2D), {len(elementos1d_msh)} Elem(1D)")

        nodes_xdmf_coords = self.mesh.geometry.x
        
        # --- MAPEAMENTO DE NÓS ---
        for msh_id, msh_coords in nodes_msh.items():
            for xdmf_id, xdmf_coords in enumerate(nodes_xdmf_coords):
                if np.linalg.norm(np.array(msh_coords) - xdmf_coords) < TOLERANCIA_COORDENADAS:
                    self.node_map_msh_to_xdmf[msh_id], self.node_map_xdmf_to_msh[xdmf_id] = xdmf_id, msh_id
                    break
        print(f"  ✅ Mapeamento de Nós criado: {len(self.node_map_xdmf_to_msh)}/{len(nodes_msh)}")

        # --- MAPEAMENTO DE ELEMENTOS 2D (CÉLULAS) ---
        cell_to_vertex = self.mesh.topology.connectivity(self.mesh.topology.dim, 0)
        # <--- CORRIGIDO: Loop seguro sobre as tags das células
        elementos2d_xdmf_centroids = []
        for idx, cell_id in enumerate(self.cell_tags.indices):
            nodes = cell_to_vertex.links(cell_id)
            centroid = np.mean(nodes_xdmf_coords[nodes], axis=0)
            pg = self.cell_tags.values[idx]
            elementos2d_xdmf_centroids.append((cell_id, pg, centroid[0], centroid[1]))
        
        for msh_id, msh_pg, node_ids in elementos2d_msh:
            if all(nid in nodes_msh for nid in node_ids):
                msh_centroid = np.mean([nodes_msh[nid] for nid in node_ids], axis=0)
                for xdmf_id, xdmf_pg, xdmf_cx, xdmf_cy in elementos2d_xdmf_centroids:
                    if msh_pg == xdmf_pg and np.sqrt((msh_centroid[0] - xdmf_cx)**2 + (msh_centroid[1] - xdmf_cy)**2) < TOLERANCIA_COORDENADAS:
                        self.elem2d_map_msh_to_xdmf[msh_id], self.elem2d_map_xdmf_to_msh[xdmf_id] = xdmf_id, msh_id
                        break
        print(f"  ✅ Mapeamento de Elementos (2D) criado: {len(self.elem2d_map_xdmf_to_msh)}/{len(elementos2d_msh)}")
        
        # --- MAPEAMENTO DE ELEMENTOS 1D (FACETAS) ---
        facet_to_vertex = self.mesh.topology.connectivity(1, 0)
        # <--- CORRIGIDO: Loop seguro sobre as tags das facetas
        elementos1d_xdmf_centroids = []
        for idx, facet_id in enumerate(self.facet_tags.indices):
            nodes = facet_to_vertex.links(facet_id)
            centroid = np.mean(nodes_xdmf_coords[nodes], axis=0)
            pg = self.facet_tags.values[idx]
            elementos1d_xdmf_centroids.append((facet_id, pg, centroid[0], centroid[1]))

        for msh_id, msh_pg, node_ids in elementos1d_msh:
            if all(nid in nodes_msh for nid in node_ids):
                msh_centroid = np.mean([nodes_msh[nid] for nid in node_ids], axis=0)
                for xdmf_id, xdmf_pg, xdmf_cx, xdmf_cy in elementos1d_xdmf_centroids:
                    if msh_pg == xdmf_pg and np.sqrt((msh_centroid[0] - xdmf_cx)**2 + (msh_centroid[1] - xdmf_cy)**2) < TOLERANCIA_COORDENADAS:
                        self.elem1d_map_msh_to_xdmf[msh_id], self.elem1d_map_xdmf_to_msh[xdmf_id] = xdmf_id, msh_id
                        break
        print(f"  ✅ Mapeamento de Facetas (1D) criado: {len(self.elem1d_map_xdmf_to_msh)}/{len(elementos1d_msh)}")


    def _mapear_nomes_para_ids(self):
        # (Esta função permanece inalterada)
        print("\nConstruindo mapeamentos dinâmicos de Nomes para IDs...")
        for camada_mat in self.config.get('camadas_material', []):
            nome = camada_mat.get('nome')
            try: id_numerico = int(nome.split('_')[-1]); self.physical_surfaces[nome] = id_numerico
            except (ValueError, IndexError): print(f"  Aviso: não foi possível extrair ID do nome '{nome}'.")
        for contorno in self.config.get('contornos', []):
            nome, id = contorno.get('nome'), contorno.get('id')
            if nome and id is not None: self.physical_lines[nome] = id
        if self.physical_surfaces: print("  🔹 Superfícies mapeadas:", dict(list(self.physical_surfaces.items())[:5]))
        if self.physical_lines: print("  🔹 Linhas mapeadas:", dict(list(self.physical_lines.items())[:5]))

    def gerar_vetor_tempo(self):
        # (Esta função permanece inalterada)
        print("\n=== 1. GERANDO VETOR DE TEMPO (LÓGICA REFINADA) ===")
        g = self.config['general']
        pts = set(np.arange(0, g['tempo_final'] + g['delta_t'], g['delta_t'])) | {0.0, g['tempo_final']}
        for c in self.config['camadas']:
            pts.add(c['birth']); pts.add(c['birth'] - g['delta_t_refinado']); pts.add(c['birth'] + g['delta_t_refinado'])
        self.vetor_tempo = np.array(sorted([p for p in pts if 0 <= p <= g['tempo_final']]))
        print(f"Vetor de tempo gerado com {len(self.vetor_tempo)} pontos.")
        return self.vetor_tempo

    def identificar_blocos_tempo(self):
        # (Esta função permanece inalterada)
        print("\n=== 2. IDENTIFICANDO BLOCOS DE TEMPO ===")
        times = sorted(set([c['birth'] for c in self.config['camadas']] + [self.config['general']['tempo_final']]))
        self.blocos_tempo = [{'id': i + 1, 'inicio': times[i], 'fim': times[i+1], 'duracao': times[i+1] - times[i]} for i in range(len(times)-1)]
        for b in self.blocos_tempo: print(f"  Bloco {b['id']}: {b['inicio']}s - {b['fim']}s ({b['duracao']/86400:.1f} dias)")
        return self.blocos_tempo

    def obter_camadas_ativas(self, tempo):
        return [c['nome'] for c in self.config['camadas'] if c['birth'] <= tempo and (c['death'] is None or tempo < c['death'])]

    def obter_physical_groups_ativos(self, camadas_ativas):
        surfaces = [self.physical_surfaces[cm['nome']] for cm in self.config.get('camadas_material',[]) if cm['camada'] in camadas_ativas and cm['nome'] in self.physical_surfaces]
        lines = [self.physical_lines[c['nome']] for c in self.config.get('contornos',[]) if c.get('nasce_com_camada') in camadas_ativas and (not c.get('desativado_pela_camada') or c.get('desativado_pela_camada') not in camadas_ativas) and c.get('nome') in self.physical_lines]
        return {'surfaces': sorted(list(set(surfaces))), 'lines': sorted(list(set(lines)))}

    def obter_elementos_e_nos(self, surface_ids=None, line_ids=None):
        # (Esta função permanece inalterada)
        res = {'elementos_dominio': [], 'nos_dominio': set(), 'elementos_contorno': [], 'nos_contorno': set()}
        if surface_ids:
            for sid in surface_ids:
                elems = self.cell_tags.find(sid)
                res['elementos_dominio'].extend(elems)
                for el in elems: res['nos_dominio'].update(self.mesh.topology.connectivity(self.mesh.topology.dim, 0).links(el))
        if line_ids:
            for lid in line_ids:
                elems = self.facet_tags.find(lid)
                res['elementos_contorno'].extend(elems)
                for el in elems: res['nos_contorno'].update(self.mesh.topology.connectivity(1, 0).links(el))
        for k in res: res[k] = sorted(list(res[k]))
        return res

    def calcular_diferencas(self, atual, anterior):
        return {'entradas': sorted(list(set(atual) - set(anterior))), 'saidas': sorted(list(set(anterior) - set(atual)))}

    def analisar_todos_blocos(self):
        # (Esta função permanece inalterada)
        print("\n=== 3. ANALISANDO TODOS OS BLOCOS DE TEMPO ===")
        for i, bloco in enumerate(self.blocos_tempo):
            print(f"\nAnalisando Bloco {bloco['id']}...")
            camadas_ativas = self.obter_camadas_ativas(bloco['inicio'])
            pgs = self.obter_physical_groups_ativos(camadas_ativas)
            elementos_nos = self.obter_elementos_e_nos(surface_ids=pgs['surfaces'], line_ids=pgs['lines'])
            diferencas = {}
            if i > 0:
                anterior = self.analise_resultados[f'bloco_{i}']
                for key in ['surfaces', 'lines']: diferencas[key] = self.calcular_diferencas(pgs[key], anterior['physical_groups'][key])
                for key in elementos_nos: diferencas[key] = self.calcular_diferencas(elementos_nos[key], anterior['elementos_nos'][key])
            self.analise_resultados[f'bloco_{i+1}'] = {'info_bloco': bloco, 'camadas_ativas': camadas_ativas, 'physical_groups': pgs, 'elementos_nos': elementos_nos, 'diferencas': diferencas}

    def gerar_relatorio(self):
        # (Esta função permanece inalterada)
        print("\n" + "="*80 + "\nRELATÓRIO COMPLETO DA ANÁLISE STAGE-WISE\n" + "="*80)
        for i, (key, data) in enumerate(self.analise_resultados.items()):
            info, pgs, en = data['info_bloco'], data['physical_groups'], data['elementos_nos']
            print(f"\nBLOCK {info['id']} ({info['inicio']}s - {info['fim']}s): Active Layers: {data['camadas_ativas']}")
            print(f"  -> Active PGs: {len(pgs['surfaces'])} surfaces, {len(pgs['lines'])} lines.")
            print(f"  -> Active Entities: {len(en['elementos_dominio'])} domain elements, {len(en['nos_dominio'])} domain nodes.")

    def converter_nos_xdmf_para_msh(self, nos_xdmf):
        return sorted([self.node_map_xdmf_to_msh.get(n) for n in nos_xdmf if n in self.node_map_xdmf_to_msh])

    def converter_elementos2d_xdmf_para_msh(self, elementos_xdmf):
        return sorted([self.elem2d_map_xdmf_to_msh.get(el) for el in elementos_xdmf if el in self.elem2d_map_xdmf_to_msh])

    def converter_elementos1d_xdmf_para_msh(self, elementos_xdmf):
        return sorted([self.elem1d_map_xdmf_to_msh.get(el) for el in elementos_xdmf if el in self.elem1d_map_xdmf_to_msh])

    def criar_resultados_versao_msh(self):
        analise_resultados_msh = {}
        for bloco_key, bloco_data in self.analise_resultados.items():
            en_xdmf = bloco_data['elementos_nos']
            en_msh = {
                'elementos_dominio': self.converter_elementos2d_xdmf_para_msh(en_xdmf['elementos_dominio']),
                'nos_dominio': self.converter_nos_xdmf_para_msh(en_xdmf['nos_dominio']),
                'elementos_contorno': self.converter_elementos1d_xdmf_para_msh(en_xdmf['elementos_contorno']),
                'nos_contorno': self.converter_nos_xdmf_para_msh(en_xdmf['nos_contorno'])
            }
            diferencas_msh = {}
            if bloco_key != 'bloco_1':
                 bloco_anterior_msh = analise_resultados_msh[f'bloco_{int(bloco_key.split("_")[-1]) - 1}']
                 for key in en_msh: diferencas_msh[key] = self.calcular_diferencas(en_msh[key], bloco_anterior_msh['elementos_nos'][key])
                 for key in ['surfaces', 'lines']: diferencas_msh[key] = bloco_data['diferencas'][key]
            
            copia_bloco = bloco_data.copy()
            copia_bloco['elementos_nos'], copia_bloco['diferencas'] = en_msh, diferencas_msh
            analise_resultados_msh[bloco_key] = copia_bloco
        return analise_resultados_msh

    def salvar_resultados_duplo(self):
        # (Esta função permanece inalterada)
        def converter(obj):
            if isinstance(obj, (np.integer, np.floating)): return obj.item()
            if isinstance(obj, np.ndarray): return obj.tolist()
            if isinstance(obj, dict): return {str(k): converter(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple, set)): return [converter(i) for i in obj]
            return obj

        info = {'yaml_file': self.yaml_file, 'xdmf_file': self.xdmf_file, 'msh_file': self.msh_file,
                'mapeamentos': {
                    'physical_surfaces': self.physical_surfaces,
                    'physical_lines': self.physical_lines
                }}
        
        resultados_xdmf = converter({'info_geral': {**info, 'tipo_numeracao': 'XDMF/H5 (para solver)'}, 'vetor_tempo': self.vetor_tempo, 'blocos_tempo': self.blocos_tempo, 'analise_resultados': self.analise_resultados})
        with open("analise_stagewise_xdmf.json", 'w', encoding='utf-8') as f: json.dump(resultados_xdmf, f, indent=2, ensure_ascii=False)
        
        analise_msh = self.criar_resultados_versao_msh()
        resultados_msh = converter({'info_geral': {**info, 'tipo_numeracao': 'MSH (para conferência)'}, 'vetor_tempo': self.vetor_tempo, 'blocos_tempo': self.blocos_tempo, 'analise_resultados': analise_msh, 'correspondencia_elementos_2d_msh_para_xdmf': self.elem2d_map_msh_to_xdmf, 'correspondencia_elementos_1d_msh_para_xdmf': self.elem1d_map_msh_to_xdmf, 'correspondencia_nos_msh_para_xdmf': self.node_map_msh_to_xdmf})
        with open("analise_stagewise_msh.json", 'w', encoding='utf-8') as f: json.dump(resultados_msh, f, indent=2, ensure_ascii=False)

        print(f"\n📄 Resultados salvos em:\n  🔧 Para o SOLVER: analise_stagewise_xdmf.json\n  ✅ Para CONFERÊNCIA: analise_stagewise_msh.json")

    def executar_analise_completa(self):
        print("INICIANDO ANÁLISE COMPLETA STAGE-WISE\n" + "="*55)
        try:
            self.gerar_vetor_tempo()
            self.identificar_blocos_tempo()
            self.analisar_todos_blocos()
            self.gerar_relatorio()
            self.salvar_resultados_duplo()
            print("\n✅ ANÁLISE CONCLUÍDA COM SUCESSO!")
        except Exception as e:
            print(f"\n❌ ERRO DURANTE A ANÁLISE: {e}")
            import traceback; traceback.print_exc()

def main():
    if len(sys.argv) != 2:
        print("❌ USO: python debug_camadas.py <nome_do_caso>\n    Exemplo: python debug_camadas.py barragem2")
        return
    
    caso = Path(sys.argv[1])
    if not caso.is_dir():
        print(f"❌ Diretório do caso não encontrado: {caso}")
        return
    
    xdmf_file, yaml_file, msh_file = caso / f"{caso.name}.xdmf", caso / f"{caso.name}.yaml", caso / f"{caso.name}.msh"
    if not all(f.exists() for f in [xdmf_file, yaml_file, msh_file]):
        print(f"❌ Arquivos essenciais não encontrados em '{caso}'. Verifique .xdmf, .yaml e .msh.")
        return
    
    analisador = AnalisadorStagewiseGenerico(str(yaml_file), str(xdmf_file), str(msh_file))
    analisador.executar_analise_completa()
    
    import os
    for tipo in ["xdmf", "msh"]:
        origem, destino = f"analise_stagewise_{tipo}.json", caso / f"{caso.name}-{tipo}.json"
        if os.path.exists(origem):
            os.replace(origem, destino); print(f"Arquivo final gerado: {destino}")

if __name__ == "__main__":
    main()