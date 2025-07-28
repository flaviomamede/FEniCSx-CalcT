# DEBUG - VARIÁVEIS DA FUNÇÃO _load_mesh()

**Data/Hora:** 2025-07-28 12:22:00
**Arquivo de malha:** barragem2/barragem2.xdmf

## 1. self.mesh

**Tipo:** <class 'dolfinx.mesh.Mesh'>
**Nome:** malha
**Dimensão:** 2
**Número de células:** 80
**Número de vértices:** 57
**Número de facetas:** 136

**Coordenadas dos vértices (primeiros 10):**
```python
[[-1.5 -3.   0. ]
 [-3.  -1.5  0. ]
 [-3.  -3.   0. ]
 [ 8.   0.   0. ]
 [ 6.5  0.   0. ]
 [ 8.  -1.5  0. ]
 [ 6.5 -1.5  0. ]
 [ 8.  -3.   0. ]
 [-1.5 -1.5  0. ]
 [ 6.5 -3.   0. ]]
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
**Número de tags:** 80
**Valores únicos:** [ 1  2  3  4  5  6  7  8  9 10]
**Valores:** [ 1  4  4  1  4  4  1  1  4  4  1  1  4  4  2  1  3  3  3  3  2  2  1  6
  3  3  2  2  6  3  3  2  2  6  6  2  5  6  6  5  5  8  6  5  5  8  6  5
  5  8  8  5  7  8  8  7  7  7 10  8  7  7 10  8  7  7 10 10  9  9 10 10
  9  9 10  9  9 10  9  9]
**Índices:** [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71
 72 73 74 75 76 77 78 79]

**Contagem por Physical Group:**
- PG 1: 8 elementos
- PG 2: 8 elementos
- PG 3: 8 elementos
- PG 4: 8 elementos
- PG 5: 8 elementos
- PG 6: 8 elementos
- PG 7: 8 elementos
- PG 8: 8 elementos
- PG 9: 8 elementos
- PG 10: 8 elementos

## 5. self.facet_tags

**Tipo:** <class 'dolfinx.mesh.MeshTags'>
**Dimensão:** 1
**Número de tags:** 40
**Valores únicos:** [11 12 13 14 15 16 17 18 19 20 21]
**Valores:** [11 11 11 11 12 11 12 11 11 11 11 12 14 12 11 11 11 13 14 20 16 13 20 16
 20 15 20 21 18 21 15 18 21 21 17 19 17 19 19 19]
**Índices:** [  1   3   4   6   7   8  11  13  18  24  26  27  30  38  40  41  51  60
  62  72  73  75  79  86  87  89  94  98  99 106 110 112 113 114 123 124
 126 131 134 135]

**Contagem por Physical Group (facetas):**
- PG 11: 12 facetas
- PG 12: 4 facetas
- PG 13: 2 facetas
- PG 14: 2 facetas
- PG 15: 2 facetas
- PG 16: 2 facetas
- PG 17: 2 facetas
- PG 18: 2 facetas
- PG 19: 4 facetas
- PG 20: 4 facetas
- PG 21: 4 facetas

## 6. LIGAÇÕES HIERÁRQUICAS

### PG-Superfícies

**PG 1:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [ 0  3  6  7 10 11 15 22]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [0, 1, 2, 8, 10, 11, 14, 15, 24]

**PG 2:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [14 20 21 26 27 31 32 35]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [10, 14, 17, 21, 22, 23, 24, 27, 28]

**PG 3:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [16 17 18 19 24 25 29 30]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [12, 13, 16, 18, 19, 20, 21, 22, 27]

**PG 4:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [ 1  2  4  5  8  9 12 13]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [3, 4, 5, 6, 7, 9, 12, 13, 16]

**PG 5:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [36 39 40 43 44 47 48 51]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [24, 27, 28, 30, 32, 33, 35, 37, 38]

**PG 6:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [23 28 33 34 37 38 42 46]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [12, 18, 25, 26, 27, 29, 31, 32, 37]

**PG 7:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [52 55 56 57 60 61 64 65]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [35, 37, 38, 40, 42, 43, 45, 46, 48]

**PG 8:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [41 45 49 50 53 54 59 63]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [29, 31, 34, 36, 37, 39, 41, 42, 46]

**PG 9:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [68 69 72 73 75 76 78 79]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [45, 46, 48, 50, 51, 53, 54, 55, 56]

**PG 10:**
  - **Elementos:** 8 elementos
  - **Índices dos elementos:** [58 62 66 67 70 71 74 77]
  - **Nós:** 9 nós únicos
  - **Índices dos nós:** [39, 41, 44, 46, 47, 49, 52, 53, 56]

### PG-Lines

**PG 11:**
  - **Facetas:** 12 facetas
  - **Índices das facetas:** [ 1  3  4  6  8 13 18 24 26 40 41 51]
  - **Nós:** 13 nós únicos
  - **Índices dos nós:** [0, 1, 2, 3, 5, 7, 9, 10, 11, 16, 17, 20, 22]

**PG 12:**
  - **Facetas:** 4 facetas
  - **Índices das facetas:** [ 7 11 27 38]
  - **Nós:** 6 nós únicos
  - **Índices dos nós:** [3, 4, 11, 12, 15, 24]

**PG 13:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [60 75]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [24, 30, 35]

**PG 14:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [30 62]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [12, 25, 29]

**PG 15:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [ 89 110]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [35, 43, 48]

**PG 16:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [73 86]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [29, 34, 39]

**PG 17:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [123 126]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [48, 50, 54]

**PG 18:**
  - **Facetas:** 2 facetas
  - **Índices das facetas:** [ 99 112]
  - **Nós:** 3 nós únicos
  - **Índices dos nós:** [39, 44, 49]

**PG 19:**
  - **Facetas:** 4 facetas
  - **Índices das facetas:** [124 131 134 135]
  - **Nós:** 5 nós únicos
  - **Índices dos nós:** [49, 52, 54, 55, 56]

**PG 20:**
  - **Facetas:** 4 facetas
  - **Índices das facetas:** [72 79 87 94]
  - **Nós:** 5 nós únicos
  - **Índices dos nós:** [29, 31, 35, 37, 38]

**PG 21:**
  - **Facetas:** 4 facetas
  - **Índices das facetas:** [ 98 106 113 114]
  - **Nós:** 5 nós únicos
  - **Índices dos nós:** [39, 41, 45, 46, 48]


## Informações Adicionais

**Comm MPI:** <mpi4py.MPI.Intracomm object at 0x72d6c5f82370>
**Rank:** 0
**Size:** 1

## Resumo

- **Malha carregada com sucesso**
- **10 Physical Groups de células**
- **11 Physical Groups de facetas**
- **Total de elementos:** 80
- **Total de vértices:** 57
