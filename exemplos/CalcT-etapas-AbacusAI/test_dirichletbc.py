#!/usr/bin/python3

"""
Teste de criau00e7u00e3o de DirichletBC com FEniCSx 0.9.0
"""

import dolfinx
import ufl
from mpi4py import MPI
import numpy as np

# Criar malha
mesh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)

# Criar espau00e7o de funu00e7u00f5es
print("Criando FunctionSpace...")
V = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))
print("FunctionSpace criado com sucesso!")

# Localizar DOFs na fronteira
print("Localizando DOFs na fronteira...")
fdim = mesh.topology.dim - 1
boundary_facets = dolfinx.mesh.locate_entities_boundary(
    mesh, fdim, lambda x: np.full(x.shape[1], True, dtype=bool)
)
boundary_dofs = dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets)
print(f"DOFs na fronteira: {len(boundary_dofs)}")

# Criar funu00e7u00e3o para valor de contorno
print("Criando funu00e7u00e3o para valor de contorno...")
u_bc = dolfinx.fem.Function(V)
u_bc.x.array[:] = 1.0

# Criar condiu00e7u00e3o de contorno
print("Criando DirichletBC...")
bc = dolfinx.fem.dirichletbc(u_bc, boundary_dofs)
print("DirichletBC criado com sucesso!")