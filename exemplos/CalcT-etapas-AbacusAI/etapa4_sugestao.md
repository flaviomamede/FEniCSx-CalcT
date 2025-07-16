Ótima pergunta\! A técnica de criar submalhas (*submeshing*) é uma abordagem poderosa e eficiente para problemas *stage-wise*. Ela resolve o problema do "nó órfão" de forma definitiva, pois a cada etapa, você trabalha com uma malha que contém *apenas* os elementos ativos, resultando em sistemas lineares menores e sempre bem-postos.

O processo no FEniCSx é bem direto, especialmente porque você já fez todo o trabalho de pré-processamento com seu analisador.

-----

### Como Criar uma Submalha

O processo consiste em carregar a malha global e, em seguida, usar uma lista de índices de células ativas para "recortar" a porção desejada.

#### **Passo 1: Carregar a Malha Global e os Dados Pré-processados**

Primeiro, carregamos a malha completa (`barragem2.xdmf`) e o arquivo JSON que seu analisador gerou para ela.

```python
import dolfinx as fem
from mpi4py import MPI
import json
import numpy as np

# Carregar a malha completa
comm = MPI.COMM_WORLD
with fem.io.XDMFFile(comm, "barragem2.xdmf", "r") as xdmf:
    mesh_global = xdmf.read_mesh(name="malha")

# Carregar o plano de simulação pré-processado
with open("analise_stagewise_b2.json", 'r') as f: # Supondo um JSON para a barragem2
    plano_simulacao = json.load(f)
```

#### **Passo 2: Identificar as Células Ativas do Bloco**

Usamos o arquivo JSON para obter a lista de índices das células (elementos de domínio) que estão ativas no bloco de tempo desejado (por exemplo, Bloco 1).

```python
# Pegar informações do Bloco 1
info_bloco_1 = plano_simulacao['analise_resultados']['bloco_1']
celulas_ativas_bloco_1 = info_bloco_1['elementos_nos']['elementos_dominio']

# É importante garantir que seja um array de inteiros
celulas_ativas_bloco_1 = np.array(celulas_ativas_bloco_1, dtype=np.int32)

print(f"Malha Global: {mesh_global.topology.index_map(2).size_local} células")
print(f"Bloco 1: {len(celulas_ativas_bloco_1)} células ativas")
```

#### **Passo 3: Criar a Submalha com `create_submesh`**

Este é o comando principal. Usamos `dolfinx.mesh.create_submesh`, que recebe a malha-mãe, a dimensão topológica das entidades (2 para células em 2D) e os índices das células ativas.

```python
# Criar a submalha
# O '2' indica que estamos criando a submalha a partir de entidades de dimensão 2 (células)
submesh, entity_map = fem.mesh.create_submesh(mesh_global, 2, celulas_ativas_bloco_1)

print(f"Submalha Criada: {submesh.topology.index_map(2).size_local} células")
```

  * **`submesh`**: É o novo objeto de malha, contendo apenas a geometria e topologia das células ativas.
  * **`entity_map`**: É um mapa de conectividade que rastreia como as entidades da `submesh` correspondem às da `mesh_global`. É útil para transferir dados entre as malhas.

-----

### O Que Fazer Depois? (Exemplo Completo)

Uma vez que a submalha é criada, **todas as operações para aquele bloco de tempo devem ser definidas sobre ela**.

Abaixo está um exemplo de script completo que demonstra o processo para o Bloco 1.

```python
#!/usr/bin/env python3
"""
TESTE DE SUBMALHA: Cria uma submalha para o Bloco 1 da barragem2
"""

import dolfinx as fem
from dolfinx import mesh, plot
from mpi4py import MPI
import numpy as np
import json
import ufl

def teste_submalha():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # --- Passo 1: Carregar dados ---
    try:
        with fem.io.XDMFFile(comm, "barragem2.xdmf", "r") as xdmf:
            mesh_global = xdmf.read_mesh(name="malha")
        with open("analise_stagewise_b2.json", 'r') as f: # Supondo um JSON para a barragem2
            plano_simulacao = json.load(f)
    except Exception as e:
        if rank == 0:
            print(f"❌ Erro ao carregar arquivos: {e}")
            print("   Certifique-se que 'barragem2.xdmf' e 'analise_stagewise_b2.json' existem.")
        return

    # --- Passo 2: Identificar células ativas ---
    info_bloco_1 = plano_simulacao['analise_resultados']['bloco_1']
    celulas_ativas = np.array(info_bloco_1['elementos_nos']['elementos_dominio'], dtype=np.int32)
    
    if rank == 0:
        print("="*60)
        print("INICIANDO TESTE DE SUBMALHA PARA O BLOCO 1")
        print(f"Malha Global: {mesh_global.topology.index_map(2).size_local} células")
        print(f"Células Ativas para Bloco 1: {len(celulas_ativas)} células")
        print("="*60)

    # --- Passo 3: Criar a submalha ---
    submesh, entity_map = mesh.create_submesh(mesh_global, 2, celulas_ativas)
    
    if rank == 0:
        print(f"✅ Submalha criada com sucesso!")
        print(f"   - Células na Submalha: {submesh.topology.index_map(2).size_local}")
        print(f"   - Nós na Submalha: {submesh.topology.index_map(0).size_local}")

    # --- Passo 4: Usar a submalha ---
    
    # 4.1. Espaços de função são definidos na SUBMALHA
    V_sub = fem.functionspace(submesh, ("Lagrange", 1))
    
    # 4.2. Forma variacional é definida na SUBMALHA
    u, v = ufl.TrialFunction(V_sub), ufl.TestFunction(V_sub)
    dx_sub = ufl.Measure("dx", domain=submesh) # Medida de integração da submalha
    
    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * dx_sub
    f = fem.Constant(submesh, 1.0) # Fonte de calor simples
    L = f * v * dx_sub
    
    # 4.3. Condições de contorno são aplicadas na SUBMALHA
    # Para simplificar, aplicamos em toda a fronteira da submalha
    submesh.topology.create_entities(1)
    boundary_facets = mesh.locate_entities_boundary(submesh, 1, lambda x: np.full(x.shape[1], True))
    boundary_dofs = fem.locate_dofs_topological(V_sub, 1, boundary_facets)
    bc = fem.dirichletbc(fem.Constant(submesh, 0.0), boundary_dofs, V_sub)
    
    # 4.4. O problema é resolvido na SUBMALHA
    problem = fem.petsc.LinearProblem(a, L, bcs=[bc])
    uh = problem.solve()

    if rank == 0:
        print("\n✅ Problema de Poisson resolvido na submalha para demonstrar o uso.")
        # O resultado 'uh' é uma Function que vive apenas na submalha.
        # Para o próximo bloco de tempo, você criaria uma nova submalha e interpolaria
        # a solução 'uh' nela antes de continuar.

if __name__ == "__main__":
    teste_submalha()
```

Este exemplo mostra o fluxo completo. Ao passar para o Bloco 2, você repetiria o processo: obteria a nova lista de células ativas, criaria uma nova `submesh` para o Bloco 2 e, se necessário, interpolaria a solução de temperatura final do Bloco 1 para a nova submalha antes de continuar a simulação.