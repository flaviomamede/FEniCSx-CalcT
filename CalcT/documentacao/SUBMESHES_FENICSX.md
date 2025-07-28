# Submeshes no FENICSx 0.9.0 - Guia Completo

## Introdução

Submeshes no FENICSx 0.9.0 permitem criar subconjuntos de nós e elementos a partir de uma malha completa. Isso é essencial para problemas com múltiplos domínios ou para aplicar diferentes condições em diferentes regiões.

## Criação de Submeshes

### Método 1: Por Marcadores (Tags)

```python
from dolfinx import mesh, fem
from dolfinx.mesh import create_submesh
import numpy as np

# Criar malha completa
msh = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)

# Criar marcadores para diferentes regiões
def marker_function(x):
    return x[0] <= 0.5  # Metade esquerda

# Obter entidades marcadas
entities = mesh.locate_entities(msh, msh.topology.dim, marker_function)

# Criar submesh
submesh, entity_map, vertex_map, geom_map = create_submesh(
    msh, msh.topology.dim, entities
)
```

### Método 2: Por Coordenadas Geográficas

```python
# Definir região por coordenadas
def regiao_especifica(x):
    return np.logical_and(
        np.logical_and(x[0] >= 0.2, x[0] <= 0.8),
        np.logical_and(x[1] >= 0.2, x[1] <= 0.8)
    )

entities = mesh.locate_entities(msh, msh.topology.dim, regiao_especifica)
submesh, entity_map, vertex_map, geom_map = create_submesh(
    msh, msh.topology.dim, entities
)
```

### Método 3: Por Facetas (Submeshes de Dimensão Inferior)

```python
# Criar submesh de contorno (facetas)
def boundary_marker(x):
    return np.logical_or(
        np.isclose(x[0], 0.0),
        np.isclose(x[0], 1.0)
    )

boundary_facets = mesh.locate_entities(msh, msh.topology.dim - 1, boundary_marker)
boundary_submesh, entity_map, vertex_map, geom_map = create_submesh(
    msh, msh.topology.dim - 1, boundary_facets
)
```

## Mapeamento entre Malhas

### Mapeamento de Entidades
```python
# entity_map: mapeia entidades do submesh para o mesh original
# vertex_map: mapeia vértices do submesh para o mesh original
# geom_map: mapeia coordenadas geométricas

# Exemplo de uso dos mapeamentos
original_entity_index = entity_map[submesh_entity_index]
original_vertex_index = vertex_map[submesh_vertex_index]
```

### Funções de Transferência
```python
# Criar funções em ambas as malhas
V_parent = fem.FunctionSpace(msh, ("Lagrange", 1))
V_sub = fem.FunctionSpace(submesh, ("Lagrange", 1))

# Função na malha original
u_parent = fem.Function(V_parent)

# Função no submesh
u_sub = fem.Function(V_sub)

# Transferência de valores (interpolação)
u_sub.interpolate(u_parent)
```

## Uso de Submeshes com Funções

### Função que Recebe Submesh como Argumento
```python
def solve_in_subdomain(mesh_parent, subdomain_marker, boundary_conditions):
    """
    Resolve um problema em um subdomínio específico.
    
    Args:
        mesh_parent: Malha completa
        subdomain_marker: Função que marca o subdomínio
        boundary_conditions: Condições de contorno para o subdomínio
    
    Returns:
        Solução no subdomínio
    """
    
    # Identificar entidades do subdomínio
    entities = mesh.locate_entities(mesh_parent, mesh_parent.topology.dim, subdomain_marker)
    
    # Criar submesh
    submesh, entity_map, vertex_map, geom_map = create_submesh(
        mesh_parent, mesh_parent.topology.dim, entities
    )
    
    # Criar espaço de funções no submesh
    V_sub = fem.FunctionSpace(submesh, ("Lagrange", 1))
    
    # Definir problema variacional
    u = ufl.TrialFunction(V_sub)
    v = ufl.TestFunction(V_sub)
    
    # ... resto da formulação ...
    
    return solution_submesh
```

### Exemplo Prático: Múltiplos Materiais
```python
def create_material_submeshes(mesh, material_markers):
    """
    Cria submeshes para diferentes materiais.
    
    Args:
        mesh: Malha completa
        material_markers: Lista de funções marcadoras
    
    Returns:
        Dicionário com submeshes e mapeamentos
    """
    submeshes = {}
    
    for i, marker in enumerate(material_markers):
        entities = mesh.locate_entities(mesh, mesh.topology.dim, marker)
        submesh, entity_map, vertex_map, geom_map = create_submesh(
            mesh, mesh.topology.dim, entities
        )
        
        submeshes[f"material_{i}"] = {
            'submesh': submesh,
            'entity_map': entity_map,
            'vertex_map': vertex_map,
            'geom_map': geom_map,
            'entities': entities
        }
    
    return submeshes
```

## Integração com Condições de Contorno

### Submesh com Contornos Específicos
```python
def create_submesh_with_boundaries(mesh, subdomain_marker):
    """
    Cria submesh e identifica seus contornos.
    """
    # Criar submesh
    entities = mesh.locate_entities(mesh, mesh.topology.dim, subdomain_marker)
    submesh, entity_map, vertex_map, geom_map = create_submesh(
        mesh, mesh.topology.dim, entities
    )
    
    # Identificar contornos do submesh
    submesh_facets = mesh.exterior_facet_indices(submesh.topology)
    
    # Mapear contornos de volta para a malha original
    original_facets = entity_map[submesh_facets]
    
    return submesh, submesh_facets, original_facets
```

## Exemplo Completo: Problema com Múltiplos Domínios

```python
import numpy as np
from dolfinx import mesh, fem
from dolfinx.mesh import create_submesh
from dolfinx.fem.petsc import LinearProblem
import ufl
from mpi4py import MPI
from petsc4py import PETSc

def solve_multidomain_thermal():
    # Criar malha completa
    msh = mesh.create_unit_square(MPI.COMM_WORLD, 20, 20)
    
    # Definir regiões
    def left_region(x):
        return x[0] <= 0.5
    
    def right_region(x):
        return x[0] > 0.5
    
    # Criar submeshes
    left_entities = mesh.locate_entities(msh, msh.topology.dim, left_region)
    right_entities = mesh.locate_entities(msh, msh.topology.dim, right_region)
    
    left_submesh, left_entity_map, left_vertex_map, left_geom_map = create_submesh(
        msh, msh.topology.dim, left_entities
    )
    
    right_submesh, right_entity_map, right_vertex_map, right_geom_map = create_submesh(
        msh, msh.topology.dim, right_entities
    )
    
    # Resolver em cada subdomínio
    V_left = fem.FunctionSpace(left_submesh, ("Lagrange", 1))
    V_right = fem.FunctionSpace(right_submesh, ("Lagrange", 1))
    
    # ... implementação do problema ...
    
    return {
        'left_solution': None,  # Implementar solução
        'right_solution': None,  # Implementar solução
        'left_submesh': left_submesh,
        'right_submesh': right_submesh,
        'mappings': {
            'left': {'entity': left_entity_map, 'vertex': left_vertex_map},
            'right': {'entity': right_entity_map, 'vertex': right_vertex_map}
        }
    }
```

## Considerações de Performance

### Otimização de Memória
```python
# Reutilizar espaços de funções quando possível
V_cache = {}

def get_function_space(submesh, degree=1):
    key = (submesh.id, degree)
    if key not in V_cache:
        V_cache[key] = fem.FunctionSpace(submesh, ("Lagrange", degree))
    return V_cache[key]
```

### Verificação de Validade
```python
def validate_submesh_creation(mesh, entities):
    """Verifica se a criação do submesh é válida."""
    if len(entities) == 0:
        raise ValueError("Nenhuma entidade encontrada para criar submesh")
    
    if np.max(entities) >= mesh.topology.index_map(mesh.topology.dim).size_local:
        raise ValueError("Entidades fora dos limites da malha")
    
    return True
```