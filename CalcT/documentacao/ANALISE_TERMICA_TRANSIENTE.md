# Análise Térmica Transiente com Geração Interna de Calor - FENICSx 0.9.0

## Introdução

A análise térmica transiente com geração interna de calor é um problema fundamental em engenharia térmica. No FENICSx 0.9.0, a formulação segue a equação do calor com termo fonte:

```
ρcₚ ∂T/∂t - ∇·(k∇T) = Q(x,t)
```

Onde:
- ρ: densidade [kg/m³]
- cₚ: capacidade térmica específica [J/(kg·K)]
- k: condutividade térmica [W/(m·K)]
- Q: geração interna de calor [W/m³]
- T: temperatura [K]
- t: tempo [s]

## Formulação Variacional

### Espaço de Funções
```python
from dolfinx import fem
from dolfinx.fem import FunctionSpace, Function

# Criar espaço de funções
V = FunctionSpace(mesh, ("Lagrange", 1))

# Funções de teste e trial
v = ufl.TestFunction(V)
T = fem.Function(V)  # Solução atual
T_old = fem.Function(V)  # Solução do passo anterior
```

### Formulação Temporal
```python
import ufl
from dolfinx.fem.petsc import LinearProblem

# Parâmetros temporais
dt = fem.Constant(mesh, PETSc.ScalarType(0.1))  # Passo de tempo

# Formulação variacional
a = rho * cp / dt * ufl.inner(T, v) * ufl.dx + k * ufl.inner(ufl.grad(T), ufl.grad(v)) * ufl.dx
L = rho * cp / dt * ufl.inner(T_old, v) * ufl.dx + Q * ufl.inner(fem.Constant(mesh, PETSc.ScalarType(1.0)), v) * ufl.dx

# Resolver
problem = LinearProblem(a, L, bcs=bcs)
T = problem.solve()
```

## Geração Interna de Calor

### Modelo Constante
```python
Q = fem.Constant(mesh, PETSc.ScalarType(1000.0))  # 1000 W/m³
```

### Modelo Espacialmente Variável
```python
# Função de geração espacialmente variável
Q_expr = fem.Expression(
    lambda x: 1000.0 * np.sin(np.pi * x[0]) * np.sin(np.pi * x[1]),
    V.element.interpolation_points()
)
Q.interpolate(Q_expr)
```

### Modelo Temporal
```python
# Geração que varia com o tempo
def update_heat_generation(t):
    Q.value = 1000.0 * np.exp(-t / 3600.0)  # Decaimento exponencial
```

## Solução Desacoplada do Tempo

### Esquema de Euler Implícito
```python
# Esquema temporal
for n in range(num_steps):
    t += dt.value
    
    # Atualizar geração de calor
    update_heat_generation(t)
    
    # Resolver sistema linear
    T = problem.solve()
    
    # Atualizar solução anterior
    T_old.x.array[:] = T.x.array[:]
```

### Esquema de Euler Explícito (para comparação)
```python
# Forma explícita (menos estável)
a_exp = rho * cp / dt * ufl.inner(T, v) * ufl.dx
L_exp = rho * cp / dt * ufl.inner(T_old, v) * ufl.dx + \
        k * ufl.inner(ufl.grad(T_old), ufl.grad(v)) * ufl.dx + \
        Q * ufl.inner(fem.Constant(mesh, PETSc.ScalarType(1.0)), v) * ufl.dx
```

## Salvamento de Resultados

### Formato XDMF para Visualização Temporal
```python
from dolfinx.io import XDMFFile

with XDMFFile(MPI.COMM_WORLD, "resultados_termicos.xdmf", "w") as file:
    file.write_mesh(mesh)
    
    for n in range(num_steps):
        # ... resolver ...
        file.write_function(T, t)
```

### Salvamento de Campos Específicos
```python
# Salvar apenas em tempos específicos
if n % 10 == 0:
    file.write_function(T, t)
```

## Exemplo Completo

```python
import numpy as np
from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
import ufl
from mpi4py import MPI
from petsc4py import PETSc

# Parâmetros do problema
L = 1.0
nx = 50
ny = 50
t_end = 10.0
dt_value = 0.1

# Criar malha
msh = mesh.create_rectangle(
    MPI.COMM_WORLD,
    [np.array([0, 0]), np.array([L, L])],
    [nx, ny],
    mesh.CellType.triangle
)

# Definir espaço de funções
V = fem.FunctionSpace(msh, ("Lagrange", 1))

# Parâmetros do material
rho = 1000.0  # kg/m³
cp = 1000.0   # J/(kg·K)
k = 50.0      # W/(m·K)

# Geração de calor
Q = fem.Constant(msh, PETSc.ScalarType(1000.0))

# Condições de contorno
T_left = fem.Constant(msh, PETSc.ScalarType(373.15))  # 100°C
T_right = fem.Constant(msh, PETSc.ScalarType(293.15))  # 20°C

# Aplicar condições de contorno
left_facets = mesh.locate_entities_boundary(msh, 1, lambda x: np.isclose(x[0], 0.0))
right_facets = mesh.locate_entities_boundary(msh, 1, lambda x: np.isclose(x[0], L))

left_dofs = fem.locate_dofs_topological(V, 1, left_facets)
right_dofs = fem.locate_dofs_topological(V, 1, right_facets)

bcs = [
    fem.dirichletbc(T_left, left_dofs),
    fem.dirichletbc(T_right, right_dofs)
]

# Funções de teste e trial
v = ufl.TestFunction(V)
T = fem.Function(V)
T_old = fem.Function(V)

# Condição inicial
T_old.interpolate(lambda x: 293.15 + 80.0 * x[0])  # Gradiente inicial

# Passo de tempo
dt = fem.Constant(msh, PETSc.ScalarType(dt_value))

# Formulação variacional
a = rho * cp / dt * ufl.inner(T, v) * ufl.dx + k * ufl.inner(ufl.grad(T), ufl.grad(v)) * ufl.dx
L = rho * cp / dt * ufl.inner(T_old, v) * ufl.dx + Q * ufl.inner(fem.Constant(msh, PETSc.ScalarType(1.0)), v) * ufl.dx

# Resolver
problem = LinearProblem(a, L, bcs=bcs)

# Loop temporal
num_steps = int(t_end / dt_value)
t = 0.0

with XDMFFile(MPI.COMM_WORLD, "resultados_termicos.xdmf", "w") as file:
    file.write_mesh(msh)
    file.write_function(T_old, t)
    
    for n in range(num_steps):
        t += dt_value
        
        # Atualizar geração de calor (exemplo: decaimento exponencial)
        Q.value = 1000.0 * np.exp(-t / 5.0)
        
        # Resolver
        T = problem.solve()
        
        # Salvar resultados
        if n % 10 == 0:
            file.write_function(T, t)
        
        # Atualizar solução anterior
        T_old.x.array[:] = T.x.array[:]
```

## Verificação de Resultados

### Conservação de Energia
```python
# Calcular energia total no sistema
energy = fem.assemble_scalar(fem.form(rho * cp * T * ufl.dx))
print(f"Energia total: {energy}")
```

### Taxa de Variação de Temperatura
```python
# Calcular taxa de variação
dT_dt = (T.x.array - T_old.x.array) / dt_value
max_dT_dt = np.max(np.abs(dT_dt))
print(f"Taxa máxima de variação: {max_dT_dt} K/s")
```