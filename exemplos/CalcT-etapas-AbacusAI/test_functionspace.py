#!/usr/bin/python3

"""
Teste de criau00e7u00e3o de FunctionSpace com FEniCSx 0.9.0
"""

import dolfinx
import ufl
from mpi4py import MPI

# Criar malha
mesh = dolfinx.mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)

# Criar espau00e7o de funu00e7u00f5es usando a funu00e7u00e3o recomendada
print("Tentando criar FunctionSpace...")
V = dolfinx.fem.functionspace(mesh, ("Lagrange", 1))
print("FunctionSpace criado com sucesso!")

# Informau00e7u00f5es sobre o espau00e7o
print(f"Dimensu00e3o do espau00e7o: {V.dofmap.index_map.size_global}")