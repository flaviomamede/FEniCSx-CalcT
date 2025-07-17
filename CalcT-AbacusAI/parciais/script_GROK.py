import dolfinx
import dolfinx.fem
import dolfinx.mesh
import dolfinx.io
import ufl
from mpi4py import MPI
import numpy as np
from petsc4py.PETSc import ScalarType

# 1. Configuração inicial
comm = MPI.COMM_WORLD

# Carregar malha completa (engloba todas as camadas)
# Suponha que a malha já foi gerada com GMSH ou outro software
# e contém marcadores para células (camadas) e facetas (contornos/interfaces)
mesh_file = "domain.xdmf"
with dolfinx.io.XDMFFile(comm, mesh_file, "r") as xdmf:
    mesh = xdmf.read_mesh(name="Grid")
    cell_tags = xdmf.read_meshtags(mesh, name="CellTags")  # Marcadores de camadas
    facet_tags = xdmf.read_meshtags(mesh, name="FacetTags")  # Marcadores de contornos

# Parâmetros do problema
num_layers = 10  # Número total de camadas
t_final = 100.0  # Tempo total (s)
num_time_steps_per_layer = 10  # Passos de tempo por camada
dt = t_final / (num_layers * num_time_steps_per_layer)  # Passo de tempo
rho = 2400.0  # Densidade (kg/m³)
cp = 1000.0  # Capacidade térmica (J/kg·K)
k = 2.0  # Condutividade térmica (W/m·K)
Q = 1000.0  # Geração de calor exotérmica (W/m³, constante por simplicidade)

# Espaço de funções
V = dolfinx.fem.FunctionSpace(mesh, ("P", 1))  # Elementos P1

# Função para marcar elementos ativos
active_cells = dolfinx.fem.FunctionSpace(mesh, ("DG", 0))  # Escalar por elemento
active_indicator = dolfinx.fem.Function(active_cells)  # 1 para ativo, 0 para inativo

# 2. Função para atualizar submalha e condições de contorno
def create_submesh_and_bcs(active_layer):
    # Marcar elementos ativos (camadas de 1 até active_layer)
    active_indicator.x.array[:] = 0.0
    for layer in range(1, active_layer + 1):
        active_cells_indices = cell_tags.indices[cell_tags.values == layer]
        active_indicator.x.array[active_cells_indices] = 1.0

    # Criar submalha com base nos elementos ativos
    active_cell_indices = np.where(active_indicator.x.array == 1.0)[0]
    submesh, entity_map, vertex_map, geom_map = dolfinx.mesh.create_submesh(
        mesh, mesh.topology.dim, active_cell_indices
    )

    # Transferir marcadores para a submalha
    submesh_cell_tags = dolfinx.mesh.meshtags(
        submesh, submesh.topology.dim, entity_map, cell_tags.values[active_cell_indices]
    )
    submesh_facet_tags = dolfinx.mesh.transfer_meshtags(
        facet_tags, submesh, entity_map, vertex_map
    )

    # Definir espaço de funções na submalha
    V_sub = dolfinx.fem.FunctionSpace(submesh, ("P", 1))

    # Condições de contorno (exemplo: Dirichlet T=0 nas faces externas)
    boundary_facets = submesh_facet_tags.indices[submesh_facet_tags.values == 100]  # Tag da face externa
    boundary_dofs = dolfinx.fem.locate_dofs_topological(V_sub, submesh.topology.dim - 1, boundary_facets)
    bc = dolfinx.fem.dirichletbc(ScalarType(0.0), boundary_dofs, V_sub)

    return submesh, V_sub, bc, submesh_cell_tags, submesh_facet_tags

# 3. Formulação variacional
def setup_variational_problem(V_sub, T_old, active_indicator_sub):
    T = ufl.TrialFunction(V_sub)
    v = ufl.TestFunction(V_sub)
    T_n = T_old  # Solução no passo anterior

    # Formas bilinear e linear
    a = (rho * cp * T * v / dt + k * ufl.inner(ufl.grad(T), ufl.grad(v))) * active_indicator_sub * ufl.dx
    L = (rho * cp * T_n * v / dt + Q * v) * active_indicator_sub * ufl.dx

    return a, L

# 4. Loop principal: blocos de tempo e camadas
T_old = dolfinx.fem.Function(V)  # Solução inicial em toda a malha
T_old.x.array[:] = 0.0  # Temperatura inicial = 0

# Salvar resultados
output_file = dolfinx.io.XDMFFile(comm, "thermal_evolution.xdmf", "w")
output_file.write_mesh(mesh)

for layer in range(1, num_layers + 1):
    print(f"Simulando camada {layer}")

    # Criar submalha e condições de contorno
    submesh, V_sub, bc, submesh_cell_tags, submesh_facet_tags = create_submesh_and_bcs(layer)

    # Transferir solução antiga para a submalha
    T_old_sub = dolfinx.fem.Function(V_sub)
    # Mapear T_old da malha completa para a submalha (interpolação)
    T_old_sub.interpolate(T_old)

    # Transferir active_indicator para a submalha
    active_indicator_sub = dolfinx.fem.Function(dolfinx.fem.FunctionSpace(submesh, ("DG", 0)))
    active_indicator_sub.x.array[:] = 1.0  # Todos os elementos da submalha são ativos

    # Configurar problema variacional
    a, L = setup_variational_problem(V_sub, T_old_sub, active_indicator_sub)

    # Montar sistema linear
    problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    T_new_sub = problem.solve()

    # Resolver para cada passo de tempo dentro da camada
    for n in range(num_time_steps_per_layer):
        t = layer * num_time_steps_per_layer * dt + n * dt
        T_new_sub = problem.solve()  # Resolver
        T_old_sub.x.array[:] = T_new_sub.x.array[:]  # Atualizar solução antiga

        # Salvar resultados na malha completa
        T_old.interpolate(T_new_sub)  # Transferir de volta para a malha completa
        output_file.write_function(T_old, t)

    # Atualizar T_old para a próxima camada
    T_old.x.array[:] = T_old.x.array[:]  # Manter solução na malha completa

output_file.close()