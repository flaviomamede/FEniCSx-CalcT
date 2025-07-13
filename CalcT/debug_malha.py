import numpy as np
from mpi4py import MPI
from dolfinx import io

def gerar_resumo_malha(xdmf_file):
    """
    Lê um arquivo de malha XDMF e imprime um resumo das tags de célula e faceta.
    """
    try:
        comm = MPI.COMM_WORLD
        
        print(f"--- Resumo da Malha para o Arquivo: {xdmf_file} ---")

        with io.XDMFFile(comm, xdmf_file, "r") as xdmf:
            # Carregar a malha e as tags
            mesh = xdmf.read_mesh(name="malha")
            # Criar entidades de dimensão 1 (facetas) necessárias para ler as tags de contorno
            mesh.topology.create_entities(1)
            cell_tags = xdmf.read_meshtags(mesh, name="malha_cells")
            facet_tags = xdmf.read_meshtags(mesh, name="malha_facets")

            # Obter os IDs únicos para domínios e contornos
            discovered_cell_tags = np.unique(cell_tags.values)
            discovered_facet_tags = np.unique(facet_tags.values)
            
            # Informações da topologia da malha
            tdim = mesh.topology.dim
            num_cells = mesh.topology.index_map(tdim).size_local
            num_facets = mesh.topology.index_map(tdim - 1).size_local
            num_nodes = mesh.topology.index_map(0).size_local


            print("\n[INFORMAÇÕES GERAIS DA MALHA]")
            print(f"Dimensão da topologia (gdim): {mesh.geometry.dim}")
            print(f"Dimensão topológica (tdim): {tdim}")
            print(f"Número total de células (elementos de domínio): {num_cells}")
            print(f"Número total de facetas (elementos de contorno): {num_facets}")
            print(f"Número total de nós: {num_nodes}")

            print("\n[PHYSICAL GROUPS - DOMÍNIOS (CELLS/SURFACES)]")
            if len(discovered_cell_tags) > 0:
                print(f"Tags de Célula (domínios) descobertas: {discovered_cell_tags}")
            else:
                print("Nenhuma tag de célula encontrada.")

            print("\n[PHYSICAL GROUPS - CONTORNOS (FACETS/LINES)]")
            if len(discovered_facet_tags) > 0:
                print(f"Tags de Faceta (contornos) descobertas: {discovered_facet_tags}")
            else:
                print("Nenhuma tag de faceta encontrada.")
        
        print("\n--- Fim do Resumo ---")

    except Exception as e:
        print(f"\nOcorreu um erro ao tentar ler a malha: {e}")
        print("Verifique se os nomes 'malha', 'malha_cells' e 'malha_facets' estão corretos.")

if __name__ == "__main__":
    # Nome do seu arquivo XDMF principal
    arquivo_malha = "barragem1.xdmf"
    gerar_resumo_malha(arquivo_malha)