# DOCUMENTAÇÃO COMPLETA - BARRAGEM FENICSx (R1)

## 📋 VISÃO GERAL

O script `barragem-Gemini-R1.py` implementa uma simulação de transferência de calor em barragens construídas em etapas (stagewise), utilizando FEniCSx para análise de elementos finitos.

### 🎯 CARACTERÍSTICAS PRINCIPAIS
- **Simulação Stage-wise**: Modela a construção da barragem em etapas temporais
- **Materiais Múltiplos**: Suporte a diferentes propriedades térmicas por região
- **Geração de Calor Exotérmica**: Modelagem de reações químicas com calor
- **Sistema Robusto**: Cascata de solvers para garantir convergência
- **Logging Completo**: Registro detalhado de execução

---

## 🏗️ ESTRUTURA HIERÁRQUICA DO CÓDIGO

### 📦 CLASSES PRINCIPAIS

#### 1. `TeeOutput` - Sistema de Logging
```
TeeOutput
├── __init__(filename)          # Inicializa sistema de log
├── write(message)              # Escreve no terminal e arquivo
├── flush()                     # Força escrita
└── close()                     # Fecha arquivo de log
```

#### 2. `SimulacaoBarragemR2` - Classe Principal
```
SimulacaoBarragem
├── __init__(config_file, json_file, log_file)
│   └── _load_config_and_analysis(config_file, json_file)    # Carrega YAML e JSON
│
├── run()                      # MÉTODO PRINCIPAL
│   ├── _setup()               # FASE DE CONFIGURAÇÃO
│   ├── _run_simulation_loop() # FASE DE SIMULAÇÃO
│   └── _finalize()            # Finalização
│
├── FASE DE CONFIGURAÇÃO (_setup)
│   ├── _load_mesh()             # Carrega malha XDMF
│   ├── _setup_materials_data()  # Carrega dados dos materiais
│   ├── _setup_function_spaces() # Cria espaços de função
│   ├── _assign_materials_to_mesh() # Atribui propriedades
│   └── _set_initial_conditions() # Aplica condições iniciais
│
├── FASE DE SIMULAÇÃO (_run_simulation_loop)
│   ├── _set_initial_temperature_for_new_elements(novos_elementos, temp_por_pg, bloco_1)
│   └── Para cada passo de tempo:
│       └── _solve_timestep(dt_val, current_time)
│           ├── [EXOTÉRMICO] _update_equivalent_time_explicitly(dt_val)
│           ├── [EXOTÉRMICO] _update_heat_generation()
│           └── _solve_temperature_equation(dt_val, current_time)
│               ├── _setup_variational_problem(dt_val, current_time)
│               ├── _get_boundary_conditions(current_time)
│               ├── _solve_with_robust_cascade(problem)
│               └── _update_state()
│       └── _save_results(...)
│
├── FUNÇÕES AUXILIARES
│   ├── _find_active_block(current_time)    # Encontra bloco ativo
│   ├── _save_results(time_step_idx, current_time) # Salva resultados
│   └── _get_fallback_bcs()                # BCs de emergência
│
└── main()                      # Ponto de entrada
```

---

## 🔄 FLUXO DE EXECUÇÃO DETALHADO

### 1️⃣ **INICIALIZAÇÃO** (`__init__`)
```python
# 1. Configura MPI e logging
# 2. Carrega arquivos de configuração
# 3. Define temperatura inicial
```

### 2️⃣ **CONFIGURAÇÃO** (`_setup`)
```python
# 1. _load_mesh() - Carrega malha XDMF
#    ├── Lê malha, células e facetas
#    └── Cria conectividade topológica
#
# 2. _setup_materials_data() - Processa YAML
#    ├── Extrai propriedades dos materiais
#    ├── Configura parâmetros exotérmicos
#    └── Define se há geração de calor
#
# 3. _setup_function_spaces() - Cria espaços
#    ├── V (CG,1) - Temperatura
#    ├── V_prop (DG,0) - Propriedades
#    └── Funções para propriedades espaciais
#
# 4. _assign_materials_to_mesh() - Atribui propriedades
#    ├── Mapeia materiais → Physical Groups
#    ├── Atribui k, ρ, cp aos elementos
#    └── Salva verificação em XDMF
#
# 5. _set_initial_conditions() - Condições iniciais
#    ├── T = T_inicial em todo domínio
#    └── teq = 0, Q_heat = 0 [EXOTÉRMICO]
```

### 3️⃣ **SIMULAÇÃO TEMPORAL** (`_run_simulation_loop`)
```python
# Para cada bloco construtivo:
#   ├── Identifica pontos temporais do bloco
#   ├── Para cada passo temporal:
#   │   ├── _solve_timestep(dt, t)
#   │   ├── _update_state()
#   │   └── _save_results() [a cada 5 passos]
#   └── Avança para próximo bloco
```

### 4️⃣ **RESOLUÇÃO DE PASSO** (`_solve_timestep`)
```python
# SE exotérmico:
#   ├── _update_equivalent_time_explicitly(dt)
#   │   ├── Calcula fator de Arrhenius
#   │   └── Atualiza teq = teq_n + dt * fator
#   │
#   └── _update_heat_generation()
#       ├── Calcula Q = f(teq, parâmetros)
#       └── Aplica modelo de Hill
#
# _solve_temperature_equation(dt, t)
#   ├── _setup_variational_problem(dt, t)
#   │   ├── Monta forma fraca
#   │   ├── Termos de massa e difusão
#   │   └── Termo fonte [EXOTÉRMICO]
#   │
#   ├── _get_boundary_conditions(t)
#   │   ├── Identifica nós inativos
#   │   ├── Aplica BCs Dirichlet
#   │   └── _get_fallback_bcs() [emergência]
#   │
#   ├── _solve_with_robust_cascade(problem)
#   │   ├── GMRES+ILU
#   │   ├── CG+HYPRE
#   │   └── LU direto
#   │
#   └── _update_state()
#       ├── T_n = T [ou Tp_n = Tp]
#       └── teq_n = teq [EXOTÉRMICO]
```

---

## 🔧 FUNÇÕES CRÍTICAS

### **Sistema de Solvers** (`_solve_with_robust_cascade`)
```python
# Cascata de solvers para garantir convergência:
# 1. GMRES + ILU (rápido, boa convergência)
# 2. CG + HYPRE (robusto, paralelo)
# 3. LU direto (garantido, mas lento)
```

### **Condições de Contorno** (`_get_boundary_conditions`)
```python
# 1. Identifica nós inativos (não construídos)
# 2. Aplica T = T_inicial nos nós inativos
# 3. Aplica BCs nos contornos físicos
# 4. Fallback: BCs em todas as facetas externas
```

### **Problema Variacional** (`_setup_variational_problem`)
```python
# Forma fraca: ∫(ρcp u v + dt θ k ∇u·∇v) dx = ∫(ρcp T_n v - dt(1-θ)k ∇T_n·∇v + dt Q v) dx
# 
# Termos:
# - Massa: ρcp u v
# - Difusão: k ∇u·∇v  
# - Fonte: Q v [EXOTÉRMICO]
# - Theta: esquema temporal (θ=0.5 = Crank-Nicolson)
```

---

## 📊 ESTRUTURA DE DADOS

### **Configuração (YAML)**
```yaml
general:
  mesh_file: "barragem1.xdmf"
  output_dir: "resultados/"
  theta: 0.5

materiais:
  - nome: "concreto"
    densidade: 2400
    condutividade_termica: 2.5
    calor_especifico: 900
    hgen:
      gera_calor: true
      par_gera_calor:
        dTadinfty: 30.0
        a_dias: 1.5
        expoente: 2.0
      EaR: 4000.0
      Tref: 20.0
```

### **Análise (JSON)**
```json
{
  "vetor_tempo": [...],
  "blocos_tempo": [...],
  "analise_resultados": {
    "bloco_1": {
      "info_bloco": {...},
      "camadas_ativas": [...],
      "physical_groups": {
        "surfaces": [...],
        "lines": [...]
      },
      "elementos_nos": {
        "elementos_dominio": [...],
        "nos_dominio": [...],
        "elementos_contorno": [...],
        "nos_contorno": [...]
      }
    }
  }
}
```

---

## 🚀 USO DO SCRIPT

### **Sintaxe**
```bash
python barragem-Gemini-R1.py <pasta_do_caso>
```

### **Estrutura de Arquivos Esperada**
```
pasta_do_caso/
├── caso.yaml          # Configuração
├── caso-xdmf.json     # Análise stage-wise
├── caso.xdmf          # Malha
├── caso.h5            # Dados da malha
└── resultados/        # Saída (criado automaticamente)
```

### **Exemplo de Execução**
```bash
python barragem-Gemini-R1.py CalcT-Gemini/barragem1
```

---

## 📈 SAÍDAS GERADAS

### **Arquivos de Resultado**
- `Temperatura_passo_XXXX.xdmf` - Campo de temperatura
- `TempoEquivalente_passo_XXXX.xdmf` - Tempo equivalente [EXOTÉRMICO]
- `GeracaoCalor_passo_XXXX.xdmf` - Geração de calor [EXOTÉRMICO]
- `propriedades_materiais.xdmf` - Verificação das propriedades

### **Log de Execução**
- `log_simulacao.md` - Registro completo da execução
- Inclui timestamps, progresso e diagnósticos

---

## ⚠️ CONSIDERAÇÕES IMPORTANTES

### **Compatibilidade FEniCSx**
- Versão 0.9.0+ requerida
- Imports corrigidos para nova API
- `functionspace` (minúsculo) em vez de `FunctionSpace`

### **Performance**
- Paralelização MPI nativa
- Cascata de solvers para robustez
- Otimização para problemas grandes

### **Limitações**
- Apenas elementos 2D suportados
- Um material por Physical Group
- Condições de contorno Dirichlet apenas

---

## 🔍 DIAGNÓSTICO DE PROBLEMAS

### **Erros Comuns**
1. **Arquivos não encontrados**: Verificar estrutura de pastas
2. **Falha de solver**: Verificar malha e parâmetros
3. **Erro de mapeamento**: Verificar Physical Groups no YAML

### **Debug**
- Log detalhado em `log_simulacao.md`
- Verificação de propriedades em XDMF
- Diagnósticos de convergência por passo

---

*Documentação gerada para `barragem-Gemini-R1.py` - Versão R1* 