#!/usr/bin/env python3
"""
Script para debugar as variáveis da função _load_mesh()
Gera arquivo mesh_readed.md com informações detalhadas das variáveis
"""

import os
import sys

# Adicionar caminhos do sistema para FEniCSx
sys.path.insert(0, '/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real/lib/python3/dist-packages')
sys.path.insert(0, '/usr/lib/python3/dist-packages')

import json
import yaml
import datetime
import numpy as np
from pathlib import Path
from dolfinx import mesh, io
from mpi4py import MPI

# Adicionar o diretório atual ao path para importar a classe
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import importlib.util
import sys

# Importar o módulo com hífens no nome
spec = importlib.util.spec_from_file_location("barragem_Gemini_R2", "barragem-Gemini-R2.py")
barragem_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(barragem_module)
SimulacaoBarragemR2 = barragem_module.SimulacaoBarragemR2

class MeshDebugger(SimulacaoBarragemR2):
    """Classe para debugar as variáveis da função _load_mesh()"""
    
    def _load_mesh(self):
        """Versão modificada de _load_mesh para debug"""
        if self.rank == 0:
            print(f"   ➡️  Carregando malha de '{self.mesh_file_path}'...")
        
        with io.XDMFFile(self.comm, str(self.mesh_file_path), "r") as xdmf:
            self.mesh = xdmf.read_mesh(name="malha")
            self.mesh.topology.create_entities(1)
            self.mesh.topology.create_connectivity(1, 2)
            self.cell_tags = xdmf.read_meshtags(self.mesh, name="malha_cells")
            self.facet_tags = xdmf.read_meshtags(self.mesh, name="malha_facets")
        
        # Gerar arquivo de debug
        self._save_mesh_debug_info()
        
        if self.rank == 0:
            print(f"   ✅ Malha '{self.mesh.name}' carregada.")
    
    def _save_mesh_debug_info(self):
        """Salva informações detalhadas das variáveis da malha"""
        if self.rank != 0:  # Apenas processo master
            return
            
        debug_file = Path(self.config_file).parent / "mesh_readed.md"
        
        with open(debug_file, 'w', encoding='utf-8') as f:
            f.write("# DEBUG - VARIÁVEIS DA FUNÇÃO _load_mesh()\n\n")
            f.write(f"**Data/Hora:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"**Arquivo de malha:** {self.mesh_file_path}\n\n")
            
            # 1. self.mesh
            f.write("## 1. self.mesh\n\n")
            f.write(f"**Tipo:** {type(self.mesh)}\n")
            f.write(f"**Nome:** {self.mesh.name}\n")
            f.write(f"**Dimensão:** {self.mesh.topology.dim}\n")
            f.write(f"**Número de células:** {self.mesh.topology.index_map(self.mesh.topology.dim).size_local}\n")
            f.write(f"**Número de vértices:** {self.mesh.topology.index_map(0).size_local}\n")
            f.write(f"**Número de facetas:** {self.mesh.topology.index_map(self.mesh.topology.dim-1).size_local}\n")
            
            # Coordenadas dos vértices
            coords = self.mesh.geometry.x
            f.write(f"\n**Coordenadas dos vértices (primeiros 10):**\n")
            f.write(f"```python\n{coords[:10]}\n```\n")
            
            # 2. self.mesh.topology.create_entities(1)
            f.write("\n## 2. self.mesh.topology.create_entities(1)\n\n")
            f.write("**Função:** Cria entidades de dimensão 1 (arestas)\n")
            f.write("**Retorno:** None (modifica a topologia internamente)\n")
            f.write("**Efeito:** Adiciona conectividade de arestas à malha\n")
            
            # 3. self.mesh.topology.create_connectivity(1, 2)
            f.write("\n## 3. self.mesh.topology.create_connectivity(1, 2)\n\n")
            f.write("**Função:** Cria conectividade entre arestas (1) e células (2)\n")
            f.write("**Retorno:** None (modifica a topologia internamente)\n")
            f.write("**Efeito:** Permite navegar de arestas para células\n")
            
            # 4. self.cell_tags
            f.write("\n## 4. self.cell_tags\n\n")
            f.write(f"**Tipo:** {type(self.cell_tags)}\n")
            f.write(f"**Dimensão:** {self.cell_tags.dim}\n")
            f.write(f"**Número de tags:** {len(self.cell_tags.indices)}\n")
            f.write(f"**Valores únicos:** {np.unique(self.cell_tags.values)}\n")
            f.write(f"**Valores:** {self.cell_tags.values}\n")
            f.write(f"**Índices:** {self.cell_tags.indices}\n")
            
            # Contagem por tag
            unique_tags, counts = np.unique(self.cell_tags.values, return_counts=True)
            f.write(f"\n**Contagem por Physical Group:**\n")
            for tag, count in zip(unique_tags, counts):
                f.write(f"- PG {tag}: {count} elementos\n")
            
            # 5. self.facet_tags
            f.write("\n## 5. self.facet_tags\n\n")
            f.write(f"**Tipo:** {type(self.facet_tags)}\n")
            f.write(f"**Dimensão:** {self.facet_tags.dim}\n")
            f.write(f"**Número de tags:** {len(self.facet_tags.indices)}\n")
            f.write(f"**Valores únicos:** {np.unique(self.facet_tags.values)}\n")
            f.write(f"**Valores:** {self.facet_tags.values}\n")
            f.write(f"**Índices:** {self.facet_tags.indices}\n")
            
            # Contagem por tag de facetas
            unique_facet_tags, facet_counts = np.unique(self.facet_tags.values, return_counts=True)
            f.write(f"\n**Contagem por Physical Group (facetas):**\n")
            for tag, count in zip(unique_facet_tags, facet_counts):
                f.write(f"- PG {tag}: {count} facetas\n")
            
            # NOVA SEÇÃO: LIGAÇÕES HIERÁRQUICAS
            f.write("\n## 6. LIGAÇÕES HIERÁRQUICAS\n\n")
            
            # DETECÇÃO AUTOMÁTICA: Separar PG por tipo de tag
            # cell_tags = PGs de células (superfícies)
            # facet_tags = PGs de facetas (linhas)
            surface_pgs = list(unique_tags)  # Todos os cell_tags são superfícies
            line_pgs = list(unique_facet_tags)  # Todos os facet_tags são linhas
            
            f.write("### PG-Superfícies (detectadas automaticamente)\n\n")
            if surface_pgs:
                for pg in sorted(surface_pgs):
                    elements = self.cell_tags.find(pg)
                    f.write(f"**PG {pg}:**\n")
                    f.write(f"  - **Elementos:** {len(elements)} elementos\n")
                    f.write(f"  - **Índices dos elementos:** {elements}\n")
                    
                    # Encontrar nós únicos dos elementos
                    unique_nodes = set()
                    for elem_idx in elements:
                        # Obter vértices do elemento
                        cell_vertices = self.mesh.topology.connectivity(self.mesh.topology.dim, 0).links(elem_idx)
                        unique_nodes.update(cell_vertices)
                    
                    f.write(f"  - **Nós:** {len(unique_nodes)} nós únicos\n")
                    f.write(f"  - **Índices dos nós:** {sorted(unique_nodes)}\n\n")
            else:
                f.write("**Nenhuma superfície detectada.**\n\n")
            
            f.write("### PG-Lines (detectadas automaticamente)\n\n")
            if line_pgs:
                for pg in sorted(line_pgs):
                    facets = self.facet_tags.find(pg)
                    f.write(f"**PG {pg}:**\n")
                    f.write(f"  - **Facetas:** {len(facets)} facetas\n")
                    f.write(f"  - **Índices das facetas:** {facets}\n")
                    
                    # Encontrar nós únicos das facetas
                    unique_nodes = set()
                    for facet_idx in facets:
                        # Obter vértices da faceta
                        facet_vertices = self.mesh.topology.connectivity(self.mesh.topology.dim-1, 0).links(facet_idx)
                        unique_nodes.update(facet_vertices)
                    
                    f.write(f"  - **Nós:** {len(unique_nodes)} nós únicos\n")
                    f.write(f"  - **Índices dos nós:** {sorted(unique_nodes)}\n\n")
            else:
                f.write("**Nenhuma linha detectada.**\n\n")
            
            # Informações adicionais
            f.write("\n## Informações Adicionais\n\n")
            f.write(f"**Comm MPI:** {self.comm}\n")
            f.write(f"**Rank:** {self.rank}\n")
            f.write(f"**Size:** {self.comm.Get_size()}\n")
            
            # Resumo
            f.write("\n## Resumo\n\n")
            f.write(f"- **Malha carregada com sucesso**\n")
            f.write(f"- **{len(unique_tags)} Physical Groups de células**\n")
            f.write(f"- **{len(unique_facet_tags)} Physical Groups de facetas**\n")
            f.write(f"- **Total de elementos:** {self.mesh.topology.index_map(self.mesh.topology.dim).size_local}\n")
            f.write(f"- **Total de vértices:** {self.mesh.topology.index_map(0).size_local}\n")
        
        print(f"   📄 Informações de debug salvas em: {debug_file}")

def main():
    if len(sys.argv) != 2:
        print(f"USO: python {sys.argv[0]} <caminho_para_pasta_do_caso>")
        sys.exit(1)
    
    pasta_caso = Path(sys.argv[1])
    caso = pasta_caso.name
    config_file = pasta_caso / f"{caso}.yaml"
    json_file = pasta_caso / f"{caso}-xdmf.json"
    
    # Criar instância do debugger
    debugger = MeshDebugger(str(config_file), str(json_file))
    
    # Executar apenas a fase de setup para carregar a malha
    try:
        debugger._setup()
        print("✅ Debug concluído com sucesso!")
    except Exception as e:
        print(f"❌ Erro durante debug: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 