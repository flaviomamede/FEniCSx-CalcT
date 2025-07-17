# Simulau00e7u00e3o Termomecau00e2nica com Construu00e7u00e3o em Camadas

Implementau00e7u00e3o Python utilizando FENICSx 0.9.0 para o problema tu00e9rmico em estruturas construu00eddas em camadas, com foco em barragens de concreto.

## Arquivos Criados

### Arquivos de Cu00f3digo

- **teste_funcionalidade.py**: Implementau00e7u00e3o de testes bau00e1sicos para verificar cada componente isolado
  - Teste bau00e1sico estacionau00e1rio
  - Teste transiente simples
  - Teste de ativau00e7u00e3o de camadas

- **exemplo_final_estavel.py**: Implementau00e7u00e3o completa e robusta
  - Classe `ConstrucaoCamadasEstavel` com todos os meu00e9todos necessau00e1rios
  - Exemplo completo de barragem com tru00eas camadas
  - Funcu00e7u00f5es auxiliares para criar malha e configurar simulau00e7u00f5es

- **testar_camadas.sh**: Script para executar todos os testes automaticamente

### Documentau00e7u00e3o

- **README.md**: Documentau00e7u00e3o geral da implementau00e7u00e3o

- **ESTRATEGIA_IMPLEMENTACAO.md**: Documentau00e7u00e3o detalhada da estratu00e9gia utilizada

## Funcionalidades Implementadas

1. **Criação e gerenciamento de malhas**
   - Tags para identificar camadas
   - Tags para condições de contorno

2. **Simulação térmica**
   - Condução de calor
   - Geração interna de calor (modelo exponencial)
   - Condições de contorno Dirichlet e Robin

3. **Construção em camadas**
   - Ativação progressiva de camadas
   - Cronograma de construção configurável

4. **Visualização**
   - Saída XDMF para visualização no ParaView
   - Arquivo de evolução temporal
   - Arquivo com campos finais

## Principais Características

- **Robustez numérica**: Formulau00e7u00e3o estável que evita problemas de convergência
- **Facilidade de uso**: Interface simples e documentada
- **Performance**: Implementau00e7u00e3o eficiente e otimizada
- **Flexibilidade**: Fácil adaptau00e7u00e3o para diferentes geometrias e cronogramas

## Dependências

- FENICSx 0.9.0
- PETSc
- NumPy
- MPI

## Execução

Para executar a implementau00e7u00e3o completa:

```bash
./testar_camadas.sh
```

## Resultados

A execuu00e7u00e3o gera arquivos XDMF que podem ser visualizados no ParaView para análise dos resultados.