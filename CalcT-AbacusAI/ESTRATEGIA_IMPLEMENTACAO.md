# Estratu00e9gia de Construu00e7u00e3o em Camadas no FENICSx

## Introduu00e7u00e3o

Este documento detalha a estratu00e9gia utilizada para implementar a simulau00e7u00e3o de construu00e7u00e3o em camadas no FENICSx. O foco estu00e1 na simulau00e7u00e3o tu00e9rmica de estruturas construu00eddas em etapas, como barragens.

## Conceitos Principais

### 1. Campo de Ativau00e7u00e3o

A abordagem central usa um **campo de ativau00e7u00e3o** que define quais partes do domu00ednio estu00e3o ativas em cada etapa da simulau00e7u00e3o. Este campo u00e9 um Function com valores binu00e1rios:

- **1.0**: Elemento ativo (participa da simulau00e7u00e3o normal)
- **0.0**: Elemento inativo (nu00e3o participa da simulau00e7u00e3o)

Este campo u00e9 utilizado para modificar as propriedades do material e os termos das equau00e7u00f5es.

### 2. Propriedades Efetivas

Ao invu00e9s de criar submeshes, usamos **propriedades efetivas** modificadas pelo campo de ativau00e7u00e3o:

```python
k_eff = ativacao * k + (1 - ativacao) * 0.001  # Condutividade mu00ednima
rho_cp_eff = ativacao * rho * cp + (1 - ativacao) * 0.001  # Capacidade mu00ednima
```

Isso faz com que elementos inativos tenham propriedades quase nulas, efetivamente removendo-os da simulau00e7u00e3o sem alterar a malha.

### 3. Gerau00e7u00e3o de Calor

Cada camada recebe um modelo de gerau00e7u00e3o de calor, geralmente uma funu00e7u00e3o exponencial:

```python
geracao_atual = Q0 * np.exp(-idade / tau) if idade > 0 else 0.0
```

Onde:
- `Q0`: Gerau00e7u00e3o inicial de calor
- `idade`: Tempo desde que a camada foi lanu00e7ada
- `tau`: Constante de tempo

### 4. Formulau00e7u00e3o Linear Estu00e1vel

A formulau00e7u00e3o utiliza `LinearProblem` ao invu00e9s de `NonlinearProblem` para garantir estabilidade:

```python
a = (rho_cp_eff / dt) * u * v * dx + k_eff * dot(grad(u), grad(v)) * dx
L = (rho_cp_eff / dt) * T_old * v * dx + geracao * v * dx
```

Condicionamento especial das propriedades garante que nu00e3o haja problemas de convergu00eancia.

## Fluxo de Simulau00e7u00e3o

1. **Inicializau00e7u00e3o**:
   - Criau00e7u00e3o da malha completa com tags para as camadas
   - Inicializau00e7u00e3o de campos (temperatura, ativau00e7u00e3o, gerau00e7u00e3o)

2. **Ativau00e7u00e3o de Camadas**:
   - Seguindo um cronograma definido, camadas su00e3o ativadas em tempos especu00edficos
   - As propriedades su00e3o ajustadas automaticamente

3. **Simulau00e7u00e3o Transiente**:
   - A cada passo de tempo, as condiu00e7u00f5es de contorno su00e3o aplicadas
   - O problema linear u00e9 resolvido para a temperatura
   - A gerau00e7u00e3o de calor u00e9 atualizada com base na idade de cada camada

## Vantagens da Abordagem

1. **Robustez**: Evita problemas de convergu00eancia
2. **Simplicidade**: Usa uma u00fanica malha durante toda a simulau00e7u00e3o
3. **Eficiu00eancia**: Evita a recriau00e7u00e3o do FunctionSpace a cada etapa
4. **Flexibilidade**: Fu00e1cil adiu00e7u00e3o de novos fenu00f4menos fu00edsicos

## Conclusu00e3o

Esta abordagem fornece uma soluu00e7u00e3o robusta para simulau00e7u00e3o de construu00e7u00e3o em camadas no FENICSx, mantendo a estabilidade numu00e9rica e permitindo a representau00e7u00e3o adequada dos fenu00f4menos fu00edsicos envolvidos.