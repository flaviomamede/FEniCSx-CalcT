# Condições Iniciais e de Contorno - Atribuição a Elementos

## Introdução

No FENICSx 0.9.0, é fundamental entender a distinção entre:
- **Domínio**: Região onde a PDE é resolvida
- **Contorno**: Fronteira do domínio onde condições de contorno são aplicadas
- **Elementos**: Entidades geométricas (triângulos, tetraedros) que compõem a malha
- **Nós**: Pontos de interpolação dos elementos

## Atribuição Robusta a Elementos vs Nós

### Problema da Interpolação

Quando atribuímos condições diretamente a elementos em vez de nós, precisamos de:
1. **Interpolação** para distribuir valores aos nós
2. **Mapeamento** entre elementos e nós
3. **Verificação** de consistência

### Método 1: Atribuição via Função de Elemento

```python
from dolfinx import fem, mesh
import numpy as np
from dolfinx.fem import Function, FunctionSpace
import ufl

# Criar malha e espaço de funções
msh = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)
V = FunctionSpace(msh, ("DG", 0))  # Elementos descontínuos para valores por elemento

# Função para atribuir valores a elementos
def assign_to_elements(element_values):
    """
    Atribui valores diretamente a elementos.
    
    Args:
        element_values: Array com valores para cada elemento
    
    Returns:
        Função definida por elemento
    """
    f = Function(V)
    f.x.array[:] = element_values
    return f

# Exemplo de uso
num_elements = msh.topology.index_map(msh.topology.dim).size_local
element_temperatures = np.full(num_elements, 293.15)  # 20°C para todos os elementos
T_element = assign_to_elements(element_temperatures)
```

### Método 2: Interpolação para Nós Contínuos

```python
# Para espaços contínuos, precisamos interpolar
def interpolate_element_to_nodes(element_function, node_space):
    """
    Interpola valores de elementos para nós.
    
    Args:
        element_function: Função definida em elementos
        node_space: Espaço de funções nos nós
    
    Returns:
        Função interpolada nos nós
    """
    node_function = Function(node_space)
    
    # Para elementos contínuos, usar interpolação
    if node_space.element.basix_element.discontinuous:
        # DG: valores diretos
        node_function.interpolate(element_function)
    else:
        # CG: interpolação ponderada
        node_function.interpolate(element_function)
    
    return node_function

# Uso
V_nodes = FunctionSpace(msh, ("Lagrange", 1))
T_nodes = interpolate_element_to_nodes(T_element, V_nodes)
```

## Condições Iniciais

### Atribuição por Elemento
```python
def set_initial_condition_by_element(mesh, temperature_function):
    """
    Define condição inicial atribuindo temperatura a elementos.
    
    Args:
        mesh: Malha do problema
        temperature_function: Função que retorna temperatura para coordenadas
    
    Returns:
        Função de condição inicial
    """
    V = FunctionSpace(mesh, ("Lagrange", 1))
    T0 = Function(V)
    
    # Obter coordenadas dos centros dos elementos
    num_elements = mesh.topology.index_map(mesh.topology.dim).size_local
    
    # Para cada elemento, calcular temperatura no centroide
    for i in range(num_elements):
        # Obter célula
        cell = mesh.geometry.dofmap.links(i)
        
        # Calcular centroide
        coords = mesh.geometry.x[cell]
        centroid = np.mean(coords, axis=0)
        
        # Atribuir temperatura
        temp = temperature_function(centroid)
        
        # Interpolação para nós do elemento
        # (implementação específica necessária)
    
    return T0

# Exemplo de uso
def initial_temp_func(coord):
    x, y = coord[0], coord[1]
    return 293.15 + 50.0 * np.sin(np.pi * x) * np.sin(np.pi * y)

T_initial = set_initial_condition_by_element(msh, initial_temp_func)
```

### Atribuição Robusta com Verificação
```python
def robust_initial_condition(mesh, element_values, node_space):
    """
    Atribuição robusta de condição inicial.
    
    Args:
        mesh: Malha do problema
        element_values: Valores por elemento
        node_space: Espaço de funções nos nós
    
    Returns:
        Função de condição inicial nos nós
    """
    
    # Verificar consistência
    num_elements = mesh.topology.index_map(mesh.topology.dim).size_local
    if len(element_values) != num_elements:
        raise ValueError(f"Número de valores ({len(element_values)}) é diferente do número de elementos ({num_elements})")
    
    # Criar função nos nós
    T0 = Function(node_space)
    
    # Para espaços contínuos, usar interpolação
    if node_space.element.basix_element.family == "Lagrange":
        # Criar função auxiliar em DG para valores por elemento
        V_DG = FunctionSpace(mesh, ("DG", 0))
        T_DG = Function(V_DG)
        T_DG.x.array[:] = element_values
        
        # Interpolar para espaço contínuo
        T0.interpolate(T_DG)
    else:
        # DG: atribuição direta
        T0.x.array[:] = element_values
    
    return T0
```

## Condições de Contorno

### Distinção entre Contornos e Domínio

```python
# Identificar contornos
def identify_boundaries(mesh):
    """
    Identifica e classifica os contornos da malha.
    
    Returns:
        Dicionário com informações sobre cada contorno
    """
    boundaries = {}
    
    # Contorno esquerdo
    left_facets = mesh.locate_entities_boundary(
        mesh, mesh.topology.dim - 1,
        lambda x: np.isclose(x[0], 0.0)
    )
    boundaries['left'] = left_facets
    
    # Contorno direito
    right_facets = mesh.locate_entities_boundary(
        mesh, mesh.topology.dim - 1,
        lambda x: np.isclose(x[0], 1.0)
    )
    boundaries['right'] = right_facets
    
    # Contorno inferior
    bottom_facets = mesh.locate_entities_boundary(
        mesh, mesh.topology.dim - 1,
        lambda x: np.isclose(x[1], 0.0)
    )
    boundaries['bottom'] = bottom_facets
    
    # Contorno superior
    top_facets = mesh.locate_entities_boundary(
        mesh, mesh.topology.dim - 1,
        lambda x: np.isclose(x[1], 1.0)
    )
    boundaries['top'] = top_facets
    
    return boundaries
```

### Atribuição de Condições de Contorno a Elementos de Contorno

```python
def assign_boundary_conditions_to_elements(mesh, boundary_values):
    """
    Atribui condições de contorno a elementos de contorno.
    
    Args:
        mesh: Malha do problema
        boundary_values: Dicionário {boundary_name: value_function}
    
    Returns:
        Lista de condições de contorno
    """
    boundaries = identify_boundaries(mesh)
    V = FunctionSpace(mesh, ("Lagrange", 1))
    bcs = []
    
    for boundary_name, value_func in boundary_values.items():
        if boundary_name in boundaries:
            facets = boundaries[boundary_name]
            dofs = fem.locate_dofs_topological(V, mesh.topology.dim - 1, facets)
            
            # Criar função de valor para o contorno
            boundary_value = Function(V)
            
            # Interpolar valores para os nós do contorno
            # (implementação específica necessária)
            
            bc = fem.dirichletbc(boundary_value, dofs)
            bcs.append(bc)
    
    return bcs
```

## Interpolação Avançada

### Interpolação entre Elementos e Nós
```python
class ElementToNodeInterpolator:
    """
    Classe para interpolação robusta entre elementos e nós.
    """
    
    def __init__(self, mesh, element_space, node_space):
        self.mesh = mesh
        self.element_space = element_space
        self.node_space = node_space
        
        # Criar mapeamentos
        self._create_mappings()
    
    def _create_mappings(self):
        """Cria mapeamentos entre elementos e nós."""
        # Implementação específica para mapeamento
        pass
    
    def element_to_nodes(self, element_function):
        """Interpola valores de elementos para nós."""
        node_function = Function(self.node_space)
        
        # Usar interpolação ponderada por área
        if self.node_space.element.basix_element.family == "Lagrange":
            # Para elementos Lagrange, usar interpolação natural
            node_function.interpolate(element_function)
        
        return node_function
    
    def nodes_to_elements(self, node_function):
        """Interpola valores de nós para elementos (média)."""
        element_function = Function(self.element_space)
        
        # Calcular média dos valores nos nós do elemento
        # Implementação específica necessária
        
        return element_function
```

## Exemplo Completo: Atribuição Robusta

```python
import numpy as np
from dolfinx import fem, mesh
from dolfinx.fem import Function, FunctionSpace, dirichletbc
import ufl
from mpi4py import MPI

def robust_thermal_setup():
    # Criar malha
    msh = mesh.create_unit_square(MPI.COMM_WORLD, 10, 10)
    
    # Espaços de funções
    V_nodes = FunctionSpace(msh, ("Lagrange", 1))
    V_elements = FunctionSpace(msh, ("DG", 0))
    
    # Condição inicial por elemento
    num_elements = msh.topology.index_map(msh.topology.dim).size_local
    initial_temps = np.linspace(293.15, 373.15, num_elements)
    
    # Interpolação para nós
    T0_nodes = robust_initial_condition(msh, initial_temps, V_nodes)
    
    # Condições de contorno
    boundary_values = {
        'left': lambda x: 373.15,
        'right': lambda x: 293.15,
        'top': lambda x: 323.15,
        'bottom': lambda x: 323.15
    }
    
    bcs = assign_boundary_conditions_to_elements(msh, boundary_values)
    
    return T0_nodes, bcs

# Uso
T_initial, boundary_conditions = robust_thermal_setup()
```

## Verificação de Consistência

```python
def verify_assignment_consistency(mesh, element_values, node_function):
    """
    Verifica consistência entre atribuição em elementos e nós.
    """
    
    # Verificar valores nos centros dos elementos
    num_elements = mesh.topology.index_map(mesh.topology.dim).size_local
    
    for i in range(num_elements):
        # Obter centroide do elemento
        cell = mesh.geometry.dofmap.links(i)
        coords = mesh.geometry.x[cell]
        centroid = np.mean(coords, axis=0)
        
        # Valor esperado
        expected = element_values[i]
        
        # Valor interpolado
        actual = node_function(centroid)
        
        # Verificar diferença
        if abs(expected - actual) > 1e-6:
            print(f"Inconsistência no elemento {i}: esperado {expected}, obtido {actual}")
    
    return True
```