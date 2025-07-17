import dolfinx
import ufl
from mpi4py import MPI
from dolfinx.mesh import create_mesh, meshtags, locate_entities, CellType
from dolfinx.fem import Function, FunctionSpace, dirichletbc, form, apply_lifting, assemble_matrix, assemble_vector, solve
from petsc4py import PETSc
import numpy as np

# 1. Crie malha global e tag por camada
mesh = load_mesh_with_all_layers()
cell_tags = meshtags(mesh, mesh.topology.dim, cell_indices, cell_values)  # Ex: 0: camada1, 1: camada2...

# 2. Parâmetros térmicos
k = 2.0     # condutividade
rho = 2400  # densidade
c = 900     # calor específico

# 3. Loop por camada (blocos de tempo)
T = Function(FunctionSpace(mesh, ("CG", 1)))  # temperatura global inicial (zero)
dt = 3600.0

for camada_id in range(num_camadas):
    # 3.1 Submesh com a nova camada ativa
    domain_cells = np.where(cell_tags.values == camada_id)[0]
    submesh, sub_to_parent = dolfinx.mesh.create_submesh(mesh, mesh.topology.dim, domain_cells)
    
    V_sub = FunctionSpace(submesh, ("CG", 1))
    
    # 3.2 Projete T anterior no novo subdomínio
    T_prev = Function(V_sub)
    T_prev.interpolate(lambda x: interpolate_from_global_T(x, T, sub_to_parent))
    
    # 3.3 Defina gerador de calor exotérmico
    Q = Function(V_sub)
    Q.interpolate(lambda x: exotermic_heat_generation(x, camada_id, t))  # função depende do tempo/camada
    
    # 3.4 Formule variáveis variacionais
    u = ufl.TrialFunction(V_sub)
    v = ufl.TestFunction(V_sub)
    
    a = (rho * c / dt) * u * v * ufl.dx + k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = (rho * c / dt) * T_prev * v * ufl.dx + Q * v * ufl.dx
    
    # 3.5 Condições de contorno atualizadas para essa camada
    bc_sub = define_dynamic_bcs(submesh, camada_id)

    # 3.6 Solução
    A = assemble_matrix(a, bcs=bc_sub)
    A.assemble()
    b = assemble_vector(L)
    apply_lifting(b, [a], [bc_sub])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    
    for bc in bc_sub:
        bc.apply(A, b)
    
    T_sub = Function(V_sub)
    solve(A, T_sub.vector, b)
    
    # 3.7 Atualiza T global (ou salva T_sub como estado da camada)
    T = update_global_temperature(T, T_sub, sub_to_parent)
