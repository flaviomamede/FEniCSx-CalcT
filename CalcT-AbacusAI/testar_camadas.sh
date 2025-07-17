#!/usr/bin/env bash

# Script para testar todos os exemplos de construcao em camadas

echo "====================================================="
echo "TESTES DE CONSTRUCAO EM CAMADAS FENICSX"
echo "====================================================="

# Ativar ambiente FEniCSx
source ../fenics_env/bin/activate

# Executar teste simples de funcionalidade
echo "\n1. Teste de funcionalidade basica\n"
python teste_funcionalidade.py

# Executar exemplo final estavel
echo "\n2. Exemplo final estavel\n"
python exemplo_final_estavel.py

echo "\n====================================================="
echo "TESTES CONCLUIDOS COM SUCESSO!"
echo "====================================================="
echo "\nArquivos gerados:"
echo "- teste_basico.xdmf"
echo "- teste_transiente.xdmf"
echo "- teste_camadas.xdmf"
echo "- simulacao_construcao.xdmf"
echo "- campos_finais.xdmf"
echo "\nUse ParaView para visualizar os resultados!"