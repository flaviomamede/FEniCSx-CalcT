from fenics import *
import matplotlib.pyplot as plt

# Defina uma malha e o espaço de funções
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Defina a solução exata para a condição de contorno
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Formule o problema variacional
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v)) * dx
L = f * v * dx

# Resolva o problema
u = Function(V)
solve(a == L, u, bc)

# Salve e visualize a solução
vtkfile = File('solution.pvd')
vtkfile << u

plot(u)
plt.show()
