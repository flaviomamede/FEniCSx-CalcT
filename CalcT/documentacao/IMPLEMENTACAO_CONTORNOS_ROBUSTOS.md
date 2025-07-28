# Implementação de Contornos Robustos com Submeshes

## Visão Geral

Esta implementação fornece um sistema robusto para aplicar condições de contorno aos elementos e interpolar para os nós, considerando cada submesh com seus respectivos physical groups definidos na malha completa.

## Arquitetura

### 1. Criação de Submeshes por Camada
- Cada camada de construção gera sua própria submesh
- Mantém mapeamentos entre malha original e submesh
- Preserva os physical groups (contornos) de cada submesh

### 2. Aplicação de Contornos
- Lê diretamente do YAML as configurações de contorno
- Aplica condições aos elementos (DG-0)
- Interpola para os nós (CG-1)
- Suporta múltiplos tipos de contorno

### 3. Sincronização
- Mantém consistência entre submeshes
- Permite ativação/desativação temporal de contornos

## Estrutura YAML para Contornos

```yaml
contornos:
  - nome: "NOME_DO_CONTORNO"
    id: 11                    # Physical Line ID do Gmsh
    nasce_com_camada: "camada_1"  # Camada onde o contorno aparece
    desativado_em: 86400     # Tempo de desativação (opcional)
    tipo: "conveccao"        # dirichlet | conveccao | fluxo | robin
    material: "concreto-ar"  # Material de interface
    h: 8.0                   # Coeficiente de transferência [W/(m²⋅K)]
    t_ext: 25.0             # Temperatura externa [°C]
    fluxo: 0.0              # Fluxo de calor [W/m²] (para tipo=fluxo)
    valor: 20.0             # Valor direto (para tipo=dirichlet)
```

## Fluxo de Execução

### 1. Inicialização
```python
# Criar submeshes para cada camada
self._create_layer_submeshes()

# Configurar espaços de funções para cada submesh
self._setup_function_spaces_with_submeshes()
```

### 2. Aplicação de Contornos por Camada
```python
# Para cada camada ativa
for layer_id in active_layers:
    bcs = self._get_boundary_conditions_for_submesh(layer_id, current_time)
    # bcs contém as condições de contorno para esta submesh
```

### 3. Resolução por Submesh
```python
# Resolver problema térmico em cada submesh
for layer_id, submesh_info in self.layer_submeshes.items():
    submesh = submesh_info['submesh']
    V_sub = self.V_layers[layer_id]
    
    # Aplicar condições de contorno específicas
    bcs = self._get_boundary_conditions_for_submesh(layer_id, current_time)
    
    # Resolver sistema linear
    problem = LinearProblem(a_sub, L_sub, bcs=bcs)
    solution = problem.solve()
```

## Tipos de Contorno Suportados

### 1. Dirichlet
```yaml
tipo: "dirichlet"
valor: 25.0  # Temperatura fixa em °C
```

### 2. Convecção
```yaml
tipo: "conveccao"
h: 8.0       # Coeficiente de convecção [W/(m²⋅K)]
t_ext: 25.0  # Temperatura externa [°C]
```

### 3. Fluxo de Calor
```yaml
tipo: "fluxo"
fluxo: 100.0  # Fluxo de calor [W/m²]
t_ambiente: 20.0  # Temperatura ambiente para referência
```

### 4. Robin (Convecção Generalizada)
```yaml
tipo: "robin"
h: 10.0      # Coeficiente de transferência
t_ext: 30.0  # Temperatura externa
```

## Mapeamento de Physical Groups

### Malha Completa → Submesh
- Cada physical line da malha original é mapeada para a submesh correspondente
- Facilita a identificação de contornos em cada camada
- Preserva a topologia dos contornos

### Exemplo de Mapeamento
```python
# Malha original tem physical lines: 11, 12, 13, 14
# Submesh da camada 1 mapeia:
#   - Physical line 11 → Contorno "ISOLAMENTO_PERFEITO"
#   - Physical line 12 → Contorno "FUNDACAO_TOPO"
#   - Physical line 13 → Contorno "FACE_MONTANTE_1"
#   - Physical line 14 → Contorno "FACE_JUSANTE_1"
```

## Ativação/Desativação Temporal

### Contornos que Nascem com a Camada
- Ativados quando a camada é construída
- Baseado no cronograma de construção

### Contornos que Podem Ser Desativados
- Desativados em tempo específico (desativado_em)
- Útil para simular remoção de formas ou mudanças de condições

## Exemplo de Uso Completo

### Configuração YAML
```yaml
# Configuração de contornos para barragem
contornos:
  # Contornos permanentes
  - nome: "BASE_ISOLADA"
    id: 11
    nasce_com_camada: "camada_1"
    tipo: "fluxo"
    h: 0.0
    t_ext: n.a

  - nome: "FACE_EXPOSTA"
    id: 12
    nasce_com_camada: "camada_1"
    tipo: "conveccao"
    h: 8.0
    t_ext: 25.0

  # Contornos temporários (formas)
  - nome: "FORMA_TEMPERATURA"
    id: 13
    nasce_com_camada: "camada_1"
    desativado_em: 172800  # 2 dias após construção
    tipo: "dirichlet"
    valor: 20.0
```

### Código de Execução
```python
# Inicializar simulador
simulator = BarragemFEniCSx(config_file="barragem1.yaml")

# Configurar submeshes
simulator._create_layer_submeshes()

# Executar simulação com contornos robustos
for current_time in time_vector:
    for layer_id in active_layers:
        bcs = simulator._get_boundary_conditions_for_submesh(layer_id, current_time)
        # Resolver para esta camada...
```

## Validação e Verificação

### Verificações Implementadas
1. **Consistência de mapeamento**: Verifica se physical groups existem
2. **Validação de valores**: Garante valores físicos válidos
3. **Verificação de cobertura**: Confere se todos os contornos são tratados
4. **Teste de interpolação**: Valida interpolação DG-0 → CG-1

### Mensagens de Debug
```python
# Mensagens informativas durante execução
print(f"   - Contorno '{contour_name}' aplicado a {len(facets)} facets")
print(f"   - Valor interpolado: {bc_value}°C")
print(f"   - Camada {layer_id}: {len(bcs)} condições de contorno aplicadas")
```

## Performance Considerations

### Otimizações
- **Cache de mapeamentos**: Mapeamentos são calculados uma vez
- **Interpolação eficiente**: Usa DG-0 → CG-1 otimizado
- **Memória compartilhada**: Reusa estruturas quando possível

### Escalabilidade
- Suporta malhas grandes (>1M elementos)
- Paralelização com MPI
- Uso de memória otimizado

## Troubleshooting

### Problemas Comuns
1. **Physical group não encontrado**: Verificar IDs no Gmsh
2. **Interpolação falha**: Verificar consistência de malhas
3. **Contorno não aplicado**: Verificar cronograma de ativação

### Ferramentas de Debug
```python
# Verificar mapeamentos
print(simulator.layer_boundary_mappings)

# Visualizar contornos
simulator._visualize_boundaries(layer_id)

# Validar condições
simulator._validate_boundary_conditions()
```