# LOG DE EXECUÇÃO - SIMULAÇÃO DE BARRAGEM (R2)
# Data: 2025-07-16 15:05:21

================================================================================
🏗️  INICIALIZANDO SIMULAÇÃO DE BARRAGEM STAGE-WISE COM FENICSx (R2)
================================================================================
🔄 Carregando configuração de 'barragem2/barragem2.yaml' e análise de 'barragem2/barragem2-xdmf.json'...
   ✅ Configuração carregada com sucesso.
   🌡️ Temperatura Inicial Padrão: 21.0°C
   📁 Diretório de Saída: barragem2/resultados
   ✅ Configuração de 3 camadas carregada.
   ✅ Configuração de 11 contornos carregada.

--- FASE DE CONFIGURAÇÃO ---
   ➡️  Carregando malha de 'barragem2/barragem2.xdmf'...
   ✅ Malha 'malha' carregada.
   ➡️  Carregando dados dos materiais do YAML...
   ➡️  Definindo espaços de função...
   ✅ Espaços de função definidos.
   ➡️  Atribuindo propriedades dos materiais aos elementos da malha...
      - Material 'fundacao' atribuído a 8 elementos no PG 1.
      - Material 'fundacao' atribuído a 8 elementos no PG 2.
      - Material 'fundacao' atribuído a 8 elementos no PG 3.
      - Material 'fundacao' atribuído a 8 elementos no PG 4.
      - Material 'concreto_face' atribuído a 8 elementos no PG 5.
      - Material 'concreto_massa' atribuído a 8 elementos no PG 6.
      - Material 'concreto_face' atribuído a 8 elementos no PG 7.
      - Material 'concreto_massa' atribuído a 8 elementos no PG 8.
      - Material 'concreto_face' atribuído a 8 elementos no PG 9.
      - Material 'concreto_massa' atribuído a 8 elementos no PG 10.
   ✅ Propriedades atribuídas. Verificação em: 'barragem2/barragem2-mat.xdmf'

--- FASE DE SIMULAÇÃO ---

❌ ERRO FATAL DURANTE A EXECUÇÃO: Function.interpolate() got an unexpected keyword argument 'cells'
