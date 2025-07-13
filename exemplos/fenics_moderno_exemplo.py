#!/usr/bin/env python3
"""
Exemplo Profissional: FEniCSx com Gmsh para Geometrias Complexas
Baseado nas melhores práticas descobertas para trabalho profissional
"""

import gmsh
import dolfinx
import dolfinx.io.gmshio
from dolfinx import fem, mesh, io
from mpi4py import MPI
import numpy as np
import ufl

def criar_geometria_complexa():
    """
    Cria geometria complexa no Gmsh com Physical Groups
    """
    # 1. Inicializar Gmsh
    gmsh.initialize()
    
    # 2. Criar modelo
    gmsh.model.add("barragem_complexa")
    
    # 3. Definir geometria (exemplo: barragem com 3 camadas)
    # Pontos base
    gmsh.model.geo.addPoint(0, 0, 0, 0.1, 1)      # Ponto 1
    gmsh.model.geo.addPoint(10, 0, 0, 0.1, 2)     # Ponto 2
    gmsh.model.geo.addPoint(10, 5, 0, 0.1, 3)     # Ponto 3
    gmsh.model.geo.addPoint(8, 7, 0, 0.1, 4)      # Ponto 4
    gmsh.model.geo.addPoint(2, 7, 0, 0.1, 5)      # Ponto 5
    gmsh.model.geo.addPoint(0, 5, 0, 0.1, 6)      # Ponto 6
    
    # Linhas
    gmsh.model.geo.addLine(1, 2, 1)  # Base
    gmsh.model.geo.addLine(2, 3, 2)  # Lateral direita
    gmsh.model.geo.addLine(3, 4, 3)  # Topo direita
    gmsh.model.geo.addLine(4, 5, 4)  # Topo
    gmsh.model.geo.addLine(5, 6, 5)  # Topo esquerda
    gmsh.model.geo.addLine(6, 1, 6)  # Lateral esquerda
    
    # Curve Loop e Surface
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4, 5, 6], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    
    # 4. Physical Groups (ESSENCIAL!)
    gmsh.model.addPhysicalGroup(1, [1], 1)  # Base - ID 1
    gmsh.model.addPhysicalGroup(1, [2, 6], 2)  # Laterais - ID 2
    gmsh.model.addPhysicalGroup(1, [3, 4, 5], 3)  # Topo - ID 3
    gmsh.model.addPhysicalGroup(2, [1], 10)  # Domínio - ID 10
    
    # Nomear Physical Groups
    gmsh.model.setPhysicalName(1, 1, "base")
    gmsh.model.setPhysicalName(1, 2, "laterais")
    gmsh.model.setPhysicalName(1, 3, "topo")
    gmsh.model.setPhysicalName(2, 10, "dominio")
    
    # 5. Sincronizar e gerar malha
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    
    # 6. Definir ordem da malha
    gmsh.model.mesh.setOrder(2)  # Malha quadrática
    
    return gmsh.model

def converter_para_dolfinx(gmsh_model, gdim=2):
    """
    Converte modelo Gmsh para DOLFINx
    """
    # Converter usando gmshio
    domain, cell_tags, facet_tags = dolfinx.io.gmshio.model_to_mesh(
        gmsh_model, MPI.COMM_WORLD, 0, gdim=gdim
    )
    
    return domain, cell_tags, facet_tags

def salvar_malha_xdmf(domain, cell_tags, facet_tags, filename="malha_complexa"):
    """
    Salva malha em formato XDMF para visualização
    """
    # Salvar malha
    with io.XDMFFile(domain.comm, f"{filename}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
    
    # Salvar tags de células
    if cell_tags is not None:
        with io.XDMFFile(domain.comm, f"{filename}_cells.xdmf", "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_meshtags(cell_tags, domain.geometry)
    
    # Salvar tags de facetas
    if facet_tags is not None:
        with io.XDMFFile(domain.comm, f"{filename}_facets.xdmf", "w") as xdmf:
            xdmf.write_mesh(domain)
            xdmf.write_meshtags(facet_tags, domain.geometry)

def problema_termico_profissional(domain, cell_tags, facet_tags):
    """
    Exemplo de problema térmico usando a malha complexa
    """
    # Espaço de elementos finitos
    from dolfinx.fem import functionspace
    V = functionspace(domain, ("CG", 2))  # Elementos quadráticos
    
    # Condições de contorno
    def base_boundary(x):
        return np.isclose(x[1], 0.0)  # y = 0
    
    def topo_boundary(x):
        return x[1] > 6.5  # y > 6.5
    
    # Marcadores de contorno
    fdim = domain.topology.dim - 1
    base_facets = mesh.locate_entities_boundary(domain, fdim, base_boundary)
    topo_facets = mesh.locate_entities_boundary(domain, fdim, topo_boundary)
    
    # Condições de contorno Dirichlet
    T_base = fem.Constant(domain, 25.0)  # 25°C na base
    T_topo = fem.Constant(domain, 20.0)  # 20°C no topo
    
    from dolfinx.fem import dirichletbc, locate_dofs_topological
    bc_base = dirichletbc(T_base, locate_dofs_topological(V, fdim, base_facets), V)
    bc_topo = dirichletbc(T_topo, locate_dofs_topological(V, fdim, topo_facets), V)
    bcs = [bc_base, bc_topo]
    
    # Variáveis do problema
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    
    # Parâmetros do material
    k = fem.Constant(domain, 2.6)  # Condutividade térmica
    Q = fem.Constant(domain, 1000.0)  # Geração de calor
    
    # Formulação fraca
    a = k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = Q * v * ufl.dx
    
    # Resolver sistema
    from dolfinx.fem.petsc import LinearProblem
    problem = LinearProblem(a, L, bcs=bcs)
    uh = problem.solve()
    
    return uh

def main():
    """
    Função principal demonstrando workflow completo
    """
    print("🚀 Iniciando workflow profissional FEniCSx + Gmsh")
    
    # 1. Criar geometria complexa
    print("📐 Criando geometria complexa...")
    gmsh_model = criar_geometria_complexa()
    
    # 2. Converter para DOLFINx
    print("🔄 Convertendo para DOLFINx...")
    domain, cell_tags, facet_tags = converter_para_dolfinx(gmsh_model)
    
    # 3. Salvar malha
    print("💾 Salvando malha em XDMF...")
    salvar_malha_xdmf(domain, cell_tags, facet_tags)
    
    # 4. Resolver problema térmico
    print("🔥 Resolvendo problema térmico...")
    solution = problema_termico_profissional(domain, cell_tags, facet_tags)
    
    # 5. Salvar resultados
    print("📊 Salvando resultados...")
    with io.XDMFFile(domain.comm, "temperatura_resultado.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(solution)
    
    # 6. Finalizar Gmsh
    gmsh.finalize()
    
    print("✅ Workflow completo!")
    print("📁 Arquivos gerados:")
    print("   - malha_complexa.xdmf (malha)")
    print("   - temperatura_resultado.xdmf (resultados)")
    print("   - Visualize no ParaView")

if __name__ == "__main__":
    main() 