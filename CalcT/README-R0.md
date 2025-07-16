Vou analisar o código `barragem-Gemini-R0.py` e listar todas as funções na ordem de execução, criando um fluxo de execução completo:

Read file: CalcT-Gemini/barragem-Gemini-R0.py
# **FLUXO DE EXECUÇÃO DO CÓDIGO `barragem-Gemini-R0.py`**

## **1. PONTO DE ENTRADA**
```
main() → SimulacaoBarragem.__init__()
```

## **2. INICIALIZAÇÃO**
```
SimulacaoBarragem.__init__() → _load_config()
```

## **3. FASE DE CONFIGURAÇÃO (método `run()`)**
```
run() → _setup() → {
    _load_mesh()
    _setup_materials()
    _setup_function_spaces()
    _set_initial_conditions()
}
```

## **4. FASE DE SIMULAÇÃO (método `run()`)**
```
run() → _run_simulation_loop() → {
    Para cada bloco de tempo:
        Para cada passo de tempo:
            _solve_timestep() → {
                _update_equivalent_time_explicitly()  # (se exotérmico)
                _update_heat_generation()            # (se exotérmico)
                _solve_temperature_equation() → {
                    _setup_variational_problem()
                    _get_boundary_conditions() → {
                        _find_active_block()
                        _get_fallback_bcs()  # (se necessário)
                    }
                    _solve_with_robust_cascade()
                }
            }
            _update_state()
            _save_results()
        }
    }
}
```

## **5. FINALIZAÇÃO**
```
run() → _finalize()
```

---

## **DETALHAMENTO DAS FUNÇÕES POR ORDEM DE EXECUÇÃO:**

### **FASE 1: INICIALIZAÇÃO**
1. `main()` - Ponto de entrada, valida argumentos
2. `SimulacaoBarragem.__init__()` - Construtor da classe
3. `_load_config()` - Carrega YAML e JSON, configura parâmetros

### **FASE 2: CONFIGURAÇÃO**
4. `_setup()` - Orquestra a configuração inicial
5. `_load_mesh()` - Carrega malha XDMF e tags
6. `_setup_materials()` - Configura propriedades dos materiais
7. `_setup_function_spaces()` - Cria espaços de função
8. `_set_initial_conditions()` - Aplica condições iniciais

### **FASE 3: SIMULAÇÃO TEMPORAL**
9. `_run_simulation_loop()` - Loop principal da simulação
10. `_solve_timestep()` - Resolve um passo de tempo
11. `_update_equivalent_time_explicitly()` - (Exotérmico) Atualiza tempo equivalente
12. `_update_heat_generation()` - (Exotérmico) Calcula geração de calor
13. `_solve_temperature_equation()` - Resolve equação da temperatura
14. `_setup_variational_problem()` - Monta problema variacional
15. `_find_active_block()` - Encontra bloco ativo no tempo atual
16. `_get_boundary_conditions()` - Define condições de contorno
17. `_solve_with_robust_cascade()` - Resolve sistema linear
18. `_update_state()` - Atualiza estado para próximo passo
19. `_save_results()` - Salva resultados periodicamente

### **FASE 4: FINALIZAÇÃO**
20. `_finalize()` - Finaliza simulação e gera relatório

---

## **FUNÇÕES AUXILIARES (chamadas conforme necessário):**
- `_get_fallback_bcs()` - BCs de emergência
- `TeeOutput.__init__()`, `write()`, `flush()`, `close()` - Logging

Este fluxo implementa uma simulação de transferência de calor stagewise com suporte a reações exotérmicas, seguindo a estrutura modular sugerida no `sugestao1.md`.