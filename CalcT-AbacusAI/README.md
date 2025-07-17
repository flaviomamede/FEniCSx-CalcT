# Construção em Camadas - FENICSx

Implementação robusta para simulação de construção em camadas no FENICSx 0.9.0

## Arquivos principais

- **teste_funcionalidade.py**: Testes básicos para verificar funcionalidades isoladas
- **exemplo_final_estavel.py**: Implementação completa e robusta de construção em camadas
- **testar_camadas.sh**: Script para executar todos os testes

## Execução

Para executar a implementação completa:

```bash
./testar_camadas.sh
```

## Características da Implementação

1. **Estratégia robusta**: Usa campos de ativação para controlar quais elementos participam da simulação
2. **Estabilidade**: Formulação que converge garantidamente (LinearProblem ao invés de NonlinearProblem)
3. **Performance**: Rápido e eficiente, com estatísticas de tempo
4. **Visualização**: Gera arquivos XDMF para visualização no ParaView

## Saídas

A execução gera os seguintes arquivos:

- `teste_basico.xdmf` - Resultado do teste básico estacionário
- `teste_transiente.xdmf` - Resultado do teste transiente
- `teste_camadas.xdmf` - Resultado do teste com ativação de camadas
- `simulacao_construcao.xdmf` - Simulação completa de barragem (evolução temporal)
- `campos_finais.xdmf` - Campos finais da simulação (temperatura, ativação, geração)

## Visualização

Os resultados podem ser visualizados utilizando o ParaView, abrindo os arquivos XDMF gerados.