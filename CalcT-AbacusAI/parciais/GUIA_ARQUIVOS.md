# Guia de Arquivos - Construção em Camadas FENICSx

## Resumo dos Arquivos Criados

### 📋 Documentação
- **`README.md`** - Guia completo de uso e implementação
- **`ESTRATEGIA_CONSTRUCAO_CAMADAS.md`** - Documentação técnica detalhada da estratégia
- **`GUIA_ARQUIVOS.md`** - Este arquivo (resumo de todos os arquivos)

### 🔧 Implementações Principais

#### 1. **`construcao_camadas_fenicsx.py`**
**Uso:** Implementação completa e robusta da classe principal
**Quando usar:** Para projetos complexos que precisam de todas as funcionalidades
**Características:**
- Classe `ConstrucaoCamadas` completa
- Suporte a submeshes dinâmicas
- Condições de contorno complexas
- Debugging avançado
- Integração com cronograma flexível

#### 2. **`simulacao_completa_camadas.py`**
**Uso:** Exemplo funcional completo pronto para executar
**Quando usar:** Para entender como tudo funciona em conjunto
**Características:**
- Classe `SimuladorConstrucaoCamadas` funcional
- Exemplo de malha com tags
- Loop de simulação completo
- Saída em XDMF
- Estatísticas de convergência

#### 3. **`exemplo_pratico_simples.py`**
**Uso:** Versão simplificada para começar rapidamente
**Quando usar:** Para primeiros testes ou implementações simples
**Características:**
- Classe `ConstrucaoCamadasSimples` minimalista
- Fácil de entender e modificar
- Poucas dependências
- Ideal para aprendizado

#### 4. **`exemplo_integrado_completo.py`**
**Uso:** Exemplo prático de barragem com todas as funcionalidades
**Quando usar:** Para ver um caso real de aplicação
**Características:**
- Classe `BarragemConstrucaoCamadas` específica
- Malha de exemplo automática
- Cronograma realista
- Relatórios detalhados
- Pronto para executar

#### 5. **`configuracao_malha_camadas.py`**
**Uso:** Geração de malhas com GMSH e configuração de tags
**Quando usar:** Para criar malhas complexas com tags apropriadas
**Características:**
- Função `criar_malha_barragem_camadas()`
- Integração com GMSH
- Tags automáticas para camadas
- Conversão para FENICSx

## Qual Arquivo Usar?

### 🚀 Para Começar Rapidamente
```python
# Use este arquivo primeiro
from exemplo_pratico_simples import ConstrucaoCamadasSimples

# Implementação mais simples
sim = ConstrucaoCamadasSimples('malha.msh', parametros)
sim.simular(cronograma)
```

### 🏗️ Para Projetos Completos
```python
# Use para implementações robustas
from construcao_camadas_fenicsx import ConstrucaoCamadas

# Implementação completa
sim = ConstrucaoCamadas(mesh, camadas_info, parametros)
sim.simular_construcao(cronograma)
```

### 📊 Para Casos Realistas
```python
# Use para ver exemplo prático
from exemplo_integrado_completo import BarragemConstrucaoCamadas

# Exemplo de barragem
barragem = BarragemConstrucaoCamadas(mesh, cell_tags, facet_tags)
barragem.simular_barragem(cronograma, parametros)
```

### 🔧 Para Criar Malhas
```python
# Use para gerar malhas com tags
from configuracao_malha_camadas import criar_malha_barragem_camadas

# Gerar malha automática
mesh, cell_tags, facet_tags = criar_malha_barragem_camadas()
```

## Fluxo de Trabalho Recomendado

### 1. **Iniciante (Primeira Vez)**
1. Leia `README.md`
2. Execute `exemplo_pratico_simples.py`
3. Modifique parâmetros para seus dados
4. Teste com sua malha

### 2. **Usuário Intermediário**
1. Leia `ESTRATEGIA_CONSTRUCAO_CAMADAS.md`
2. Execute `exemplo_integrado_completo.py`
3. Use `configuracao_malha_camadas.py` para sua geometria
4. Adapte `simulacao_completa_camadas.py` para seu caso

### 3. **Usuário Avançado**
1. Use `construcao_camadas_fenicsx.py` como base
2. Implemente funcionalidades específicas
3. Integre com análise termomecânica
4. Otimize para casos grandes

## Dependências por Arquivo

### Mínimas (exemplo_pratico_simples.py)
```python
import numpy as np
import dolfinx
import dolfinx.fem as fem
import dolfinx.io
import dolfinx.nls.petsc
import ufl
from mpi4py import MPI
```

### Completas (demais arquivos)
```python
# Adicione às dependências mínimas:
from petsc4py import PETSc
import time
import gmsh  # Para configuracao_malha_camadas.py
```

## Estrutura de Dados

### Parâmetros Termicos
```python
parametros = {
    'k': 2.5,      # Condutividade [W/m·K]
    'rho': 2400.0, # Densidade [kg/m³]
    'cp': 1000.0   # Calor específico [J/kg·K]
}
```

### Configuração de Camadas
```python
camadas_info = [
    {
        'id': 0,
        'tempo_lancamento': 0.0,
        'material': {'Q0': 1000.0, 'tau': 48.0},
        'dirichlet': [{'tag': 1, 'valor': 20.0}],
        'robin': [{'tag': 2, 'h': 10.0, 'T_amb': 15.0}]
    }
]
```

### Cronograma
```python
cronograma = [
    {'tempo': 0.0, 'tipo': 'nova_camada', 'camada': 0},
    {'tempo': 168.0, 'tipo': 'nova_camada', 'camada': 1},
    {'tempo': 336.0, 'tipo': 'fim'}
]
```

## Arquivos de Saída

### Gerados Automaticamente
- **`resultado.xdmf`** - Evolução temporal da temperatura
- **`campos_finais.xdmf`** - Campos finais para análise
- **`debug.xdmf`** - Campos de ativação e geração para debug

### Visualização
Use ParaView ou similar para visualizar arquivos XDMF:
```bash
paraview resultado.xdmf
```

## Customização

### Para Sua Geometria
1. Modifique `configuracao_malha_camadas.py`
2. Ajuste tags das camadas
3. Defina condições de contorno apropriadas

### Para Seus Materiais
1. Modifique `parametros_termicos`
2. Ajuste leis de geração de calor
3. Inclua propriedades específicas

### Para Seu Cronograma
1. Defina tempos de lançamento
2. Especifique condições por camada
3. Configure critérios de parada

## Solução de Problemas

### Problemas Comuns
1. **Malha sem tags**: Use `configuracao_malha_camadas.py`
2. **Convergência lenta**: Reduza `dt` ou ajuste solver
3. **Temperaturas irreais**: Verifique geração de calor
4. **Erro de tags**: Verifique numeração das tags (1, 2, 3, ...)

### Debug
```python
# Salvar campos para debug
sim.salvar_campos_debug("debug.xdmf")

# Verificar ativação
print(f"Camadas ativas: {sim.camadas_ativas}")
print(f"DOFs ativos: {np.sum(sim.ativacao.x.array > 0.5)}")
```

## Próximos Passos

### Extensões Possíveis
1. **Acoplamento termomecânico**
2. **Múltiplos materiais**
3. **Condições climáticas**
4. **Otimização de performance**

### Melhorias Sugeridas
1. Interface gráfica para configuração
2. Pós-processamento automático
3. Análise de sensibilidade
4. Validação experimental

---

**Escolha o arquivo mais adequado ao seu nível de experiência e necessidades específicas!**