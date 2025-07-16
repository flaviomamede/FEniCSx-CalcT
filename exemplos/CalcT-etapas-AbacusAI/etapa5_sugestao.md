A etapa 4 funcionou. Conseguimos rodar sobre uma submesh.


Agora a etapa 5 é rodar uma submesh após a outra.


Temos duas abordagens possíveis. Ou:


1) a submesh seguinte é uma submesh que engloba a anterior mais os elementos que nascem

2) a submesh seguinte é só submesh dos elementos que nascem


Em ambos os casos o difícil é fazer com que as variáveis de estado (teq, Q e T) da submesh anterior continuem de onde pararam, sendo que os novos elementos nascidos recebem as condições iniciais. 

RESPOSTA:

Excelente\! O sucesso da Etapa 4 é o marco mais importante que alcançamos. Conseguir rodar a simulação em uma única submalha significa que todos os componentes individuais — a física, a formulação variacional e o solver — estão corretos.

Você está certíssimo ao identificar a Etapa 5 e seu principal desafio. Sua análise das duas abordagens é perspicaz, e a dificuldade que você apontou é exatamente o "coração" da simulação *stage-wise* com submalhas.

Vamos estruturar a solução definitiva.

-----

### A Abordagem Correta: Submalha Cumulativa

Entre as duas opções que você propôs, a **Abordagem 1 (Submalha Cumulativa)** é o caminho correto e fisicamente representativo para um problema de condução de calor. A cada nova etapa construtiva, o novo bloco de concreto troca calor com todos os blocos já existentes.

A Abordagem 2 (submalha apenas com elementos novos) não funcionaria, pois trataria cada nova camada como um corpo termicamente isolado, o que não reflete a realidade da condução de calor na barragem.

Portanto, nosso foco será implementar a **submalha cumulativa**.

-----

### O Desafio: Transferência de Estado entre Malhas

Seu diagnóstico está perfeito: o problema é transferir as variáveis de estado (`T`, `teq`) de uma submalha para a próxima. Como a malha muda a cada bloco, os Espaços de Função também mudam. Você não pode simplesmente copiar os arrays de valores, pois eles têm tamanhos e mapeamentos de graus de liberdade (DOFs) diferentes.

A ferramenta correta no FEniCSx para resolver isso é a **interpolação**.

O método `Function.interpolate()` permite "projetar" uma função de uma malha de origem para uma malha de destino. Ele avalia a função de origem nos pontos da malha de destino para calcular os novos valores dos coeficientes.

-----

### Plano de Jogo para a Etapa 5: A Simulação Stage-Wise Completa

Aqui está o pseudo-código e a lógica para construir o solver final.

```plaintext
# =================================================================
# ALGORITMO PARA A ETAPA 5 - SOLVER STAGE-WISE COMPLETO
# =================================================================

1.  INICIALIZAÇÃO GERAL:
        Carregar malha global e dados do JSON (como na Etapa 4).
        Inicializar variável para guardar a solução final de cada bloco:
        `solucao_T_bloco_anterior = None`
        `solucao_teq_bloco_anterior = None`

2.  LOOP EXTERNO (Iterar sobre os blocos construtivos):
        PARA CADA `bloco` EM `plano_simulacao['blocos_tempo']`:

        2.1. CRIAR A SUBMALHA CUMULATIVA ATUAL:
            Obter a lista de TODAS as células ativas até o bloco atual.
            (Ex: para o Bloco 2, a lista será `celulas_bloco1 + celulas_bloco2`).
            Chamar `submesh, entity_map = fem.mesh.create_submesh(...)` com essa lista cumulativa.

        2.2. DEFINIR ESPAÇOS E FUNÇÕES NA NOVA SUBMALHA:
            Definir `V_T`, `V_teq`, etc., na `submesh` recém-criada.
            Criar as `Function`s para este bloco: `T_atual`, `T_anterior`, `teq_atual`, `teq_anterior`.

        2.3. TRANSFERIR ESTADO (O PASSO CRÍTICO):
            SE `solucao_T_bloco_anterior` NÃO for Nula (ou seja, não estamos no primeiro bloco):
                // Interpolar a solução final do bloco anterior para a condição inicial deste bloco.
                `T_anterior.interpolate(solucao_T_bloco_anterior)`
                `teq_anterior.interpolate(solucao_teq_bloco_anterior)`
            SENÃO (estamos no primeiro bloco):
                // Usar a condição inicial global.
                `T_anterior.x.array[:] = T_inicial`
                `teq_anterior.x.array[:] = 0.0`

        2.4. EXECUTAR LOOP TEMPORAL INTERNO (Como na Etapa 4):
            Pegar o vetor de tempo correspondente APENAS a este bloco.
            PARA CADA `passo_de_tempo` DENTRO DO BLOCO ATUAL:
                - Calcular `dt`.
                - `update_equivalent_time()` para obter `teq_atual`.
                - `calculate_heat_generation()` para obter `Q_atual`.
                - Montar e resolver o `LinearProblem` para obter `T_atual`.
                - Atualizar `T_anterior.x.array[:] = T_atual.x.array[:]`.
                - Atualizar `teq_anterior.x.array[:] = teq_atual.x.array[:]`.
                - Salvar resultados (opcional).

        2.5. PREPARAR PARA O PRÓXIMO BLOCO:
            Salvar as funções finais deste bloco para serem usadas na próxima iteração do loop externo.
            `solucao_T_bloco_anterior = T_atual`
            `solucao_teq_bloco_anterior = teq_atual`

3.  FIM DO LOOP EXTERNO. SIMULAÇÃO COMPLETA.
```

### Esqueleto de Código Python

Para te dar um ponto de partida, aqui está um esqueleto de como o código principal poderia ser estruturado.

```python
import dolfinx as fem
from mpi4py import MPI
# ... outros imports

def run_full_stagewise_simulation():
    # --- PASSO 1: INICIALIZAÇÃO GERAL ---
    comm = MPI.COMM_WORLD
    mesh_global, plano_simulacao = carregar_dados_iniciais() # Sua função para carregar

    solucao_T_bloco_anterior = None
    solucao_teq_bloco_anterior = None

    todos_elementos_ativos = []

    # --- PASSO 2: LOOP EXTERNO SOBRE OS BLOCOS ---
    for bloco_key, info_bloco in plano_simulacao['analise_resultados'].items():
        if comm.rank == 0:
            print(f"\n🚀 INICIANDO BLOCO CONSTRUTIVO: {bloco_key} 🚀")

        # --- 2.1: Criar Submalha Cumulativa ---
        novos_elementos = info_bloco['elementos_nos']['elementos_dominio']
        todos_elementos_ativos.extend(novos_elementos)
        
        celulas_cumulativas = np.array(sorted(list(set(todos_elementos_ativos))), dtype=np.int32)
        
        submesh, entity_map = fem.mesh.create_submesh(mesh_global, 2, celulas_cumulativas)

        # --- 2.2: Definir Espaços e Funções na Nova Submalha ---
        V_T = fem.functionspace(submesh, ("Lagrange", 1))
        V_teq = fem.functionspace(submesh, ("Lagrange", 1))
        # ... criar T_atual, T_anterior, teq_atual, teq_anterior na submalha ...
        T_anterior = fem.Function(V_T)
        teq_anterior = fem.Function(V_teq)
        
        # --- 2.3: Transferir Estado Anterior ---
        if solucao_T_bloco_anterior is not None:
            if comm.rank == 0:
                print(f"   Interpolando estado do bloco anterior...")
            T_anterior.interpolate(solucao_T_bloco_anterior)
            teq_anterior.interpolate(solucao_teq_bloco_anterior)
            
            # TODO: Aplicar condição inicial nos nós novos, se necessário.
            # Uma forma é identificar os nós da submalha que não existiam antes.
            
        else: # Primeiro bloco
            if comm.rank == 0:
                print(f"   Aplicando condição inicial global...")
            T_anterior.x.array[:] = 20.0  # T_inicial
            teq_anterior.x.array[:] = 0.0

        # --- 2.4: Executar Loop Temporal Interno ---
        # (Aqui entra a lógica do seu main() da Etapa 4, adaptada para
        # usar as funções e espaços da submalha atual)
        T_final_deste_bloco, teq_final_deste_bloco = executar_loop_interno(submesh, T_anterior, teq_anterior, info_bloco)
        
        # --- 2.5: Salvar Estado para Próximo Bloco ---
        solucao_T_bloco_anterior = T_final_deste_bloco
        solucao_teq_bloco_anterior = teq_final_deste_bloco

    if comm.rank == 0:
        print("\n🎉 SIMULAÇÃO STAGE-WISE COMPLETA CONCLUÍDA! 🎉")

# Funções auxiliares a serem implementadas
def carregar_dados_iniciais():
    # ... seu código para carregar malha e json ...
    pass

def executar_loop_interno(submesh, T_n, teq_n, info_bloco):
    # ... lógica de loop da Etapa 4 (cálculo de dt, teq, Q, solve) ...
    # ... retorna as funções de T e teq no final do bloco ...
    pass
```

Este plano de ação e o esqueleto de código te dão um caminho claro para construir a Etapa 5. O segredo está no loop externo que gerencia as submalhas e na interpolação para garantir a continuidade entre as etapas.