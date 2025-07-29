# DEBUG - VARIÁVEIS DA FUNÇÃO _load_mesh()

**Data/Hora:** 2025-07-28 18:35:47
**Arquivo de malha:** barragem1/barragem1.xdmf

## 1. self.mesh

**Tipo:** <class 'dolfinx.mesh.Mesh'>
**Nome:** malha
**Dimensão:** 2
**Número de células:** 10
**Número de vértices:** 19
**Número de facetas:** 28

**Coordenadas dos vértices (primeiros 10):**
```python
[[ 0.  0.  0.]
 [-3.  0.  0.]
 [ 0. -3.  0.]
 [-3. -3.  0.]
 [ 1.  0.  0.]
 [ 1. -3.  0.]
 [ 8.  0.  0.]
 [ 5.  0.  0.]
 [ 8. -3.  0.]
 [ 5. -3.  0.]]
```

## 2. self.mesh.topology.create_entities(1)

**Função:** Cria entidades de dimensão 1 (arestas)
**Retorno:** None (modifica a topologia internamente)
**Efeito:** Adiciona conectividade de arestas à malha

## 3. self.mesh.topology.create_connectivity(1, 2)

**Função:** Cria conectividade entre arestas (1) e células (2)
**Retorno:** None (modifica a topologia internamente)
**Efeito:** Permite navegar de arestas para células

## 4. self.cell_tags

**Tipo:** <class 'dolfinx.mesh.MeshTags'>
**Dimensão:** 2
**Número de tags:** 10
**Valores únicos:** [ 1  2  3  4  5  6  7  8  9 10]
**Valores:** [ 1  2  4  5  3  6  7  8  9 10]
**Índices:** [0 1 2 3 4 5 6 7 8 9]

**Contagem por Physical Group:**
- PG 1: 1 elementos
- PG 2: 1 elementos
- PG 3: 1 elementos
- PG 4: 1 elementos
- PG 5: 1 elementos
- PG 6: 1 elementos
- PG 7: 1 elementos
- PG 8: 1 elementos
- PG 9: 1 elementos
- PG 10: 1 elementos

## 5. self.facet_tags

**Tipo:** <class 'dolfinx.mesh.MeshTags'>
**Dimensão:** 1
**Número de tags:** 20
**Valores únicos:** [11 12 13 14 15 16 17 18 19 20 21]
**Valores:** [12 13 11 11 11 11 12 11 14 11 20 20 15 16 21 21 17 18 19 19]
**Índices:** [ 0  3  4  5  6 10 11 12 14 15 16 17 19 20 21 22 24 25 26 27]

**Contagem por Physical Group (facetas):**
- PG 11: 6 facetas
- PG 12: 2 facetas
- PG 13: 1 facetas
- PG 14: 1 facetas
- PG 15: 1 facetas
- PG 16: 1 facetas
- PG 17: 1 facetas
- PG 18: 1 facetas
- PG 19: 2 facetas
- PG 20: 2 facetas
- PG 21: 2 facetas

## 6. LIGAÇÕES HIERÁRQUICAS

### PG-Superfícies (detectadas automaticamente)

**PG 1:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [0]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [0, 1, 2, 3]

**PG 2:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [1]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [0, 2, 4, 5]

**PG 3:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [4]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [4, 5, 7, 9]

**PG 4:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [2]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [6, 7, 8, 9]

**PG 5:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [3]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [0, 4, 10, 11]

**PG 6:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [5]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [4, 7, 10, 12]

**PG 7:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [6]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [10, 11, 13, 14]

**PG 8:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [7]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [10, 12, 13, 15]

**PG 9:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [8]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [13, 14, 16, 17]

**PG 10:**
  - **Elementos:** 1 elementos
  - **Índices dos elementos:** [9]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [13, 15, 16, 18]

### PG-Lines (detectadas automaticamente)

**PG 11:**
  - **Facetas:** 6 facetas
  - **Índices das facetas:** [ 4  5  6 10 12 15]
  - **Nós:** 7 nós únicos
  - **Índices dos nós:** [1, 2, 3, 5, 6, 8, 9]

**PG 12:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [ 0 11]
  - **Nós:** 4 nós únicos
  - **Índices dos nós:** [0, 1, 6, 7]

**PG 13:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [3]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [0, 11]

**PG 14:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [14]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [7, 12]

**PG 15:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [19]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [11, 14]

**PG 16:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [20]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [12, 15]

**PG 17:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [24]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [14, 17]

**PG 18:**
  - **Facetas:** 1 facetas
  - **Índices das facetas:** [25]
  - **Nós:** 2 nós únicos
  - **Índices dos nós:** [15, 18]

**PG 19:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [26 27]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [16, 17, 18]

**PG 20:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [16 17]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [10, 11, 12]

**PG 21:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [21 22]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [13, 14, 15]


## Informações Adicionais

**Comm MPI:** <mpi4py.MPI.Intracomm object at 0x75638b18a4c0>
**Rank:** 0
**Size:** 1

## Resumo

- **Malha carregada com sucesso**
- **10 Physical Groups de células**
- **11 Physical Groups de facetas**
- **Total de elementos:** 10
- **Total de vértices:** 19
