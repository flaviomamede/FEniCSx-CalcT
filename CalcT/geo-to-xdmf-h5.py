import gmsh
from dolfinx.io import gmshio, XDMFFile
from mpi4py import MPI
import os
import sys

def convert_geo_to_xdmf(geo_file_path, output_mesh_name="malha_barragem"):
    """
    Converte um arquivo .geo do Gmsh para o formato .xdmf para uso com FEniCSx.

    Args:
        geo_file_path (str): Caminho para o arquivo .geo do Gmsh.
        output_mesh_name (str): Nome base para os arquivos de saída (.xdmf e .h5).
    """
    # Inicializa o Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1) # Exibir saída do Gmsh no terminal

    try:
        # Cria um novo modelo Gmsh
        model_name = os.path.splitext(os.path.basename(geo_file_path))[0]
        gmsh.model.add(model_name)

        # Abre e mescla o arquivo .geo
        gmsh.open(geo_file_path)

        # Gera a malha 2D
        gmsh.model.mesh.generate(2) # 2 para malha 2D

        # Define a ordem da malha (opcional, mas recomendado para maior precisão)
        # gmsh.model.mesh.setOrder(1) # Ordem 1 para elementos lineares

        # Converte o modelo Gmsh para uma malha DOLFINx
        # gdim=2 para uma malha 2D, pois o Gmsh sempre gera pontos 3D
        mesh, cell_tags, facet_tags = gmshio.model_to_mesh(
            gmsh.model, MPI.COMM_WORLD, 0, gdim=2
        )
        mesh.name = "malha"
        cell_tags.name = f"{mesh.name}_cells"
        facet_tags.name = f"{mesh.name}_facets"

        # Salva a malha e as tags em arquivos XDMF
        output_xdmf_path = f"{output_mesh_name}.xdmf"
        output_h5_path = f"{output_mesh_name}.h5"

        with XDMFFile(MPI.COMM_WORLD, output_xdmf_path, "w") as xdmf:
            xdmf.write_mesh(mesh)
            xdmf.write_meshtags(cell_tags, mesh.geometry)
            xdmf.write_meshtags(facet_tags, mesh.geometry)

        print(f"Malha e tags salvas em {output_xdmf_path} e {output_h5_path}")

    except Exception as e:
        print(f"Ocorreu um erro durante a conversão: {e}")
    finally:
        # Finaliza o Gmsh
        gmsh.finalize()

def main():
    """
    Função principal que processa argumentos de linha de comando
    """
    if len(sys.argv) < 2:
        print("❌ USO: python geo-to-xdmf-h5.py <arquivo.geo> [nome_base_saida]")
        print("    Exemplo: python geo-to-xdmf-h5.py barragem1.geo barragem1")
        return
    
    geo_file = sys.argv[1]
    
    # Nome base para saída (sem extensão)
    if len(sys.argv) >= 3:
        output_base = sys.argv[2]
        if output_base.endswith('.xdmf'):
            output_base = output_base.replace('.xdmf', '')
        if output_base.endswith('.h5'):
            output_base = output_base.replace('.h5', '')
    else:
        # Usar nome do arquivo .geo como base
        output_base = os.path.splitext(geo_file)[0]
    
    # Verificar se arquivo existe
    if not os.path.exists(geo_file):
        print(f"❌ Arquivo não encontrado: {geo_file}")
        return
    
    print(f"🔧 Convertendo: {geo_file} → {output_base}.xdmf/.h5")
    
    # Chama a função para converter o arquivo
    convert_geo_to_xdmf(geo_file, output_base)

if __name__ == "__main__":
    main()