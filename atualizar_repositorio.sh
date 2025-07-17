#!/usr/bin/env bash

# Script para atualizar o repositório Git com as novas implementações

echo "=============================================================="
echo "ATUALIZANDO REPOSITÓRIO GIT COM NOVOS ARQUIVOS"
echo "=============================================================="

# Adicionar arquivos modificados no diretório CalcT
git add CalcT/README-R1.md
git add CalcT/barragem-Gemini-R1.py
git add CalcT/barragem2/barragem2.yaml
git add CalcT/barragem2/log_simulacao.md

# Adicionar novos arquivos no diretório CalcT
git add CalcT/README-R2.md
git add CalcT/barragem-Gemini-R2.py
git add CalcT/barragem2/log_simulacao_R2.md
git add CalcT/sugestao2.md

# Adicionar todo o novo diretório CalcT-AbacusAI
git add CalcT-AbacusAI/

# Remover arquivo deletado
git rm RESUMO_MIGRACAO_FENICSX.md

# Fazer o commit das alterações
git commit -m "Adiciona implementação robusta de construção em camadas usando FENICSx"

# Explicação do que foi feito
echo ""
echo "Commit preparado com as seguintes alterações:"
echo "- Adicionado diretório CalcT-AbacusAI/ com implementação completa"
echo "- Atualizados arquivos em CalcT/ com melhorias"
echo "- Removido arquivo RESUMO_MIGRACAO_FENICSX.md"
echo ""

# Mensagem para fazer push
echo "Para enviar as alterações para o GitHub, execute:"
echo "git push origin main"
echo ""
echo "=============================================================="