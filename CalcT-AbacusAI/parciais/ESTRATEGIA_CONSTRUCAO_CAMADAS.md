# Estratégia para Construção em Camadas no FENICSx

## Visão Geral

Este documento descreve uma estratégia robusta para implementar problemas de construção em camadas com ativação/desativação de elementos no FENICSx, especificamente para análise térmica com geração interna de calor exotérmica (típica de barragens de concreto).

## Conceitos Fundamentais

### 1. Malha Global Única com Tags

**Estratégia**: Usar uma malha global que contém todas as camadas desde o início, mas controlada por tags.

```python
# Cada camada tem uma tag específica
camada_0 = tag 1
camada_1 = tag 2  
camada_2 = tag 3
```

**Vantagens**:
- Evita remalhamento durante a simulação
- Mantém conectividade entre camadas
- Facilita condições de contorno nas interfaces

### 2. Campos de Ativação

**Conceito**: Usar campos escalares que controlam quando elementos estão "ativos":

```python
# Campo binário: 1 = ativo, 0 = inativo
campo_ativacao = fem.Function(V)

# Propriedades efetivas
k_efetiva = campo_ativacao * k_material
rho_cp_efetiva = campo_ativacao * rho_cp_material
```

### 3. Submeshes Dinâmicas

**Aplicação**: Embora se mantenha a malha global, as propriedades são aplicadas apenas onde necessário através dos campos de ativação.

## Implementação Detalhada

### Estrutura da Classe Principal

```python
class SimuladorConstrucaoCamadas:
    def __init__(self, mesh, cell_tags, facet_tags, camadas_info, parametros):
        # Malha global
        self.mesh = mesh
        self.cell_tags = cell_tags  # Tags das células por camada
        self.facet_tags = facet_tags  # Tags das facetas para contorno
        
        # Campos de solução
        self.T = fem.Function(V)      # Temperatura atual
        self.T_old = fem.Function(V)  # Temperatura anterior
        
        # Campos de controle
        self.ativacao_total = fem.Function(V)  # Ativação global
        self.geracao_total = fem.Function(V)   # Geração de calor
        
        # Estado da simulação
        self.camadas_ativas = set()
        self.tempo_atual = 0.0
```

### Ativação de Camadas

```python
def ativar_camada(self, num_camada):
    """
    Ativa uma nova camada no tempo correto
    """
    # Adicionar à lista de camadas ativas
    self.camadas_ativas.add(num_camada)
    
    # Atualizar campos de ativação
    self.atualizar_campos_ativacao()
    
    # Invalidar solver (condições de contorno mudaram)
    self.solver = None
```

### Formulação Variacional

```python
def criar_formulacao_variacional(self, dt):
    v = ufl.TestFunction(self.V)
    
    # Propriedades efetivas com ativação
    k_eff = self.ativacao_total * k_base
    rho_cp_eff = self.ativacao_total * rho_cp_base
    
    # Formulação
    F = (rho_cp_eff * (self.T - self.T_old) / dt) * v * self.dx
    F += k_eff * ufl.dot(ufl.grad(self.T), ufl.grad(v)) * self.dx
    F -= self.geracao_total * v * self.dx
    
    # Condições de contorno dinâmicas
    for camada in self.camadas_ativas:
        # Aplicar condições específicas da camada
        ...
    
    return F
```

### Geração de Calor Exotérmica

```python
def atualizar_geracao_calor(self):
    """
    Atualiza geração baseada na idade de cada camada
    """
    for num_camada in self.camadas_ativas:
        camada_info = self.camadas_info[num_camada]
        
        # Calcular idade da camada
        idade = self.tempo_atual - camada_info['tempo_lancamento']
        
        # Lei exponencial para geração
        Q0 = camada_info['material']['Q0']
        tau = camada_info['material']['tau']
        geracao = Q0 * np.exp(-idade / tau)
        
        # Aplicar nas células da camada
        cells_camada = self.obter_cells_camada(num_camada)
        for cell in cells_camada:
            dofs = self.V.dofmap.cell_dofs(cell)
            self.geracao_total.x.array[dofs] = geracao
```

## Aspectos Críticos

### 1. Condições de Contorno Dinâmicas

```python
def definir_condicoes_dirichlet(self):
    bcs = []
    
    for num_camada in self.camadas_ativas:
        camada_info = self.camadas_info[num_camada]
        
        # Condições específicas da camada
        if 'dirichlet' in camada_info:
            for bc_info in camada_info['dirichlet']:
                facets = self.facet_tags.find(bc_info['tag'])
                dofs = fem.locate_dofs_topological(self.V, dim-1, facets)
                bc = fem.dirichletbc(bc_info['valor'], dofs, self.V)
                bcs.append(bc)
    
    return bcs
```

### 2. Interfaces Entre Camadas

- **Automaticamente tratadas**: Quando camadas adjacentes são ativadas, a continuidade é mantida automaticamente
- **Sem necessidade de condições especiais**: A formulação de Galerkin garante continuidade naturalmente

### 3. Solver Reconfiguração

```python
def resolver_passo_tempo(self, dt):
    # Atualizar campos (geração varia com tempo)
    self.atualizar_campos_ativacao()
    
    # Reconfigurar solver se nova camada foi ativada
    if self.solver is None:
        self.configurar_solver(dt)
    
    # Resolver sistema
    n_iter, converged = self.solver.solve(self.T)
    
    return converged
```

## Cronograma de Construção

```python
cronograma = [
    {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},
    {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},  # 7 dias
    {'tempo': 336.0, 'tipo': 'nova_camada', 'camada': 2},  # 14 dias
    {'tempo': 720.0, 'tipo': 'fim_simulacao'}               # 30 dias
]
```

## Integração com GMSH

### Criação de Malha com Tags

```python
def criar_malha_barragem():
    gmsh.initialize()
    
    # Criar geometria por camadas
    for i in range(num_camadas):
        # Criar superfície da camada
        surface = gmsh.model.geo.addPlaneSurface([loop])
        
        # Adicionar tag física
        gmsh.model.addPhysicalGroup(2, [surface], i+1)
        gmsh.model.setPhysicalName(2, i+1, f"camada_{i}")
    
    # Gerar malha
    gmsh.model.mesh.generate(2)
    
    # Converter para FENICSx
    mesh, cell_tags, facet_tags = dolfinx.io.gmshio.read_from_msh(
        "barragem.msh", MPI.COMM_WORLD
    )
    
    return mesh, cell_tags, facet_tags
```

## Vantagens da Estratégia

1. **Robustez**: Malha única evita problemas de remalhamento
2. **Flexibilidade**: Fácil adição/remoção de camadas
3. **Eficiência**: Solver reconfigurado apenas quando necessário
4. **Generalidade**: Funciona para 2D e 3D
5. **Compatibilidade**: Totalmente compatível com FENICSx 0.9.0

## Considerações de Performance

### Otimizações Recomendadas

1. **Campos de ativação**: Usar `np.int32` para economizar memória
2. **Solver linear**: Usar MUMPS ou PARDISO para sistemas grandes
3. **Adaptação temporal**: Reduzir `dt` quando nova camada é ativada
4. **Salvamento**: Salvar apenas campos necessários para reduzir I/O

### Escalabilidade

```python
# Para problemas grandes, considerar:
opts["pc_type"] = "gamg"          # Multigrid algébrico
opts["ksp_type"] = "cg"           # Gradiente conjugado
opts["pc_gamg_agg_nsmooths"] = 1  # Otimização GAMG
```

## Debugging e Diagnóstico

### Campos de Debug

```python
def salvar_campos_debug(self):
    # Salvar campos de ativação
    self.ativacao_total.name = "Ativacao"
    self.geracao_total.name = "Geracao_Calor"
    
    arquivo.write_function(self.ativacao_total, tempo)
    arquivo.write_function(self.geracao_total, tempo)
```

### Verificações de Consistência

```python
def verificar_consistencia(self):
    # Verificar se todas as camadas ativas têm células
    for camada in self.camadas_ativas:
        cells = self.obter_cells_camada(camada)
        if len(cells) == 0:
            print(f"Aviso: Camada {camada} não tem células!")
    
    # Verificar continuidade de temperatura
    T_min, T_max = np.min(self.T.x.array), np.max(self.T.x.array)
    if T_max > 200.0:  # Temperatura muito alta
        print(f"Aviso: Temperatura máxima = {T_max:.1f}°C")
```

## Extensões Futuras

### 1. Análise Termomecânica

```python
# Adicionar deslocamentos
V_mech = fem.VectorFunctionSpace(mesh, ("Lagrange", 1))
u = fem.Function(V_mech)

# Acoplamento termo-mecânico
epsilon_th = alpha * (T - T_ref) * ufl.Identity(dim)
```

### 2. Múltiplos Materiais

```python
# Diferentes materiais por camada
def obter_propriedades_camada(self, num_camada):
    material = self.camadas_info[num_camada]['material']
    return material['k'], material['rho'], material['cp']
```

### 3. Condições de Contorno Evoluindo

```python
# Condições que mudam com o tempo
def atualizar_condicoes_contorno(self, tempo):
    # Exemplo: convecção que muda com clima
    h_convectivo = self.calcular_conveccao_climatica(tempo)
    return h_convectivo
```

## Resumo

Esta estratégia oferece uma solução robusta e eficiente para problemas de construção em camadas no FENICSx, mantendo a flexibilidade necessária para análises complexas de barragens de concreto com geração interna de calor.