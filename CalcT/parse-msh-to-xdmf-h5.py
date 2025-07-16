import gmsh
from dolfinx.io import gmshio, XDMFFile
from mpi4py import MPI
import os
import sys

def convert_msh_to_xdmf(msh_file_path, output_mesh_name="malha_barragem"):
    """
    Converte um arquivo .msh do Gmsh para o formato .xdmf para uso com FEniCSx.

    Args:
        msh_file_path (str): Caminho para o arquivo .msh do Gmsh.
        output_mesh_name (str): Nome base para os arquivos de saída (.xdmf e .h5).
    """
    # Inicializa o Gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1) # Exibir saída do Gmsh no terminal

    try:
        # Abre diretamente o arquivo .msh
        gmsh.open(msh_file_path)

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
        print("❌ USO: python parse-msh-to-xdmf-h5.py <diretorio_caso>")
        print("    Exemplo: python parse-msh-to-xdmf-h5.py barragem1")
        return
    
    # Novo: argumento é o diretório do caso
    case_dir = sys.argv[1]
    if not os.path.isdir(case_dir):
        print(f"❌ Diretório não encontrado: {case_dir}")
        return

    # Nome base = nome do diretório
    base_name = os.path.basename(os.path.normpath(case_dir))
    msh_file = os.path.join(case_dir, f"{base_name}.msh")
    output_base = os.path.join(case_dir, base_name)

    # Verificar se arquivo .msh existe
    if not os.path.exists(msh_file):
        print(f"❌ Arquivo .msh não encontrado: {msh_file}")
        return
    
    print(f"🔧 Convertendo: {msh_file} → {output_base}.xdmf/.h5")
    
    # Chama a função para converter o arquivo .msh
    convert_msh_to_xdmf(msh_file, output_base)

if __name__ == "__main__":
    main()