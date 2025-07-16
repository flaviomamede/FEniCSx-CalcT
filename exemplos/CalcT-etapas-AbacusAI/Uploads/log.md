# LOG DE EXECUÇÃO - 2025-07-13 18:17:46

🏗️  BARRAGEM STAGE-WISE - VERSÃO FINAL OTIMIZADA
============================================================
🌡️  Temperatura inicial: 20.0°C
📊 Blocos de tempo: 3
🎯 EXECUTANDO APENAS PRIMEIRO BLOCO DE TEMPO
📁 Malha: barragem1.xdmf
   ✅ Facet tags carregadas do XDMF: [11 12 13 14 15 16 17 18 19 20 21]
📁 Malha carregada: barragem1.xdmf
🏷️  Surfaces: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
🏷️  Lines: [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
✅ Materiais configurados: ['concreto_massa', 'concreto_face', 'fundacao']
🔥 Sistema exotérmico detectado!
✅ Espaços de função configurados (sistema exotérmico: Tp + teq + Q)
   ✅ Condições iniciais exotérmicas: Tp=20.0°C, teq=0s, Q=0

🚀 SIMULAÇÃO APENAS NO PRIMEIRO BLOCO DE TEMPO
   📅 Bloco 1: 0.0h → 48.0h
   🕒 Duração: 48.0h
   📊 Pontos de tempo: 23
============================================================

📅 Passo 1: t=0.0h → 0.2h (dt=0.2h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [600.0, 600.0]s
   ✅ Geração de calor calculada: [1.4e+03, 1.4e+03] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 2: t=0.2h → 2.4h (dt=2.2h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [8640.0, 8640.0]s
   ✅ Geração de calor calculada: [2.9e+05, 2.9e+05] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 3: t=2.4h → 4.8h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [17280.0, 17280.0]s
   ✅ Geração de calor calculada: [1.1e+06, 1.1e+06] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 4: t=4.8h → 7.2h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [25920.0, 25920.0]s
   ✅ Geração de calor calculada: [2.5e+06, 2.5e+06] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 5: t=7.2h → 9.6h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [34560.0, 34560.0]s
   ✅ Geração de calor calculada: [4.3e+06, 4.3e+06] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)
   💾 Salvo: temperatura_bloco1_0005.xdmf + tempo_equiv_bloco1_0005.xdmf

📅 Passo 6: t=9.6h → 12.0h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [43200.0, 43200.0]s
   ✅ Geração de calor calculada: [6.5e+06, 6.5e+06] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 7: t=12.0h → 14.4h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [51840.0, 51840.0]s
   ✅ Geração de calor calculada: [8.9e+06, 8.9e+06] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 8: t=14.4h → 16.8h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [60480.0, 60480.0]s
   ✅ Geração de calor calculada: [1.2e+07, 1.2e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 9: t=16.8h → 19.2h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [69120.0, 69120.0]s
   ✅ Geração de calor calculada: [1.4e+07, 1.4e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 10: t=19.2h → 21.6h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [77760.0, 77760.0]s
   ✅ Geração de calor calculada: [1.7e+07, 1.7e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)
   💾 Salvo: temperatura_bloco1_0010.xdmf + tempo_equiv_bloco1_0010.xdmf

📅 Passo 11: t=21.6h → 24.0h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [86400.0, 86400.0]s
   ✅ Geração de calor calculada: [2.0e+07, 2.0e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 12: t=24.0h → 26.4h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [95040.0, 95040.0]s
   ✅ Geração de calor calculada: [2.3e+07, 2.3e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 13: t=26.4h → 28.8h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [103680.0, 103680.0]s
   ✅ Geração de calor calculada: [2.5e+07, 2.5e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 14: t=28.8h → 31.2h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [112320.0, 112320.0]s
   ✅ Geração de calor calculada: [2.8e+07, 2.8e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 15: t=31.2h → 33.6h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [120960.0, 120960.0]s
   ✅ Geração de calor calculada: [3.0e+07, 3.0e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)
   💾 Salvo: temperatura_bloco1_0015.xdmf + tempo_equiv_bloco1_0015.xdmf

📅 Passo 16: t=33.6h → 36.0h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [129600.0, 129600.0]s
   ✅ Geração de calor calculada: [3.2e+07, 3.2e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 17: t=36.0h → 38.4h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [138240.0, 138240.0]s
   ✅ Geração de calor calculada: [3.4e+07, 3.4e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 18: t=38.4h → 40.8h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [146880.0, 146880.0]s
   ✅ Geração de calor calculada: [3.6e+07, 3.6e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 19: t=40.8h → 43.2h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [155520.0, 155520.0]s
   ✅ Geração de calor calculada: [3.8e+07, 3.8e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 20: t=43.2h → 45.6h (dt=2.4h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [164160.0, 164160.0]s
   ✅ Geração de calor calculada: [4.0e+07, 4.0e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)
   💾 Salvo: temperatura_bloco1_0020.xdmf + tempo_equiv_bloco1_0020.xdmf

📅 Passo 21: t=45.6h → 47.8h (dt=2.2h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [172200.0, 172200.0]s
   ✅ Geração de calor calculada: [4.1e+07, 4.1e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

📅 Passo 22: t=47.8h → 48.0h (dt=0.2h)
   🔄 Calculando tempo equivalente para 6 domínios...
   ✅ Tempo equivalente atualizado: [172800.0, 172800.0]s
   ✅ Geração de calor calculada: [4.1e+07, 4.1e+07] W/m³
   🔥 Formulação exotérmica com Q para 6 domínios ativos
   ✅ Formulação exotérmica com Q construída para 6 domínios
   🎯 CORREÇÃO: Nós ativos: 13 de 19
   🎯 CORREÇÃO: Contornos ativos: [11, 12, 13, 14, 20]
   🔧 CORREÇÃO: Nós inativos (órfãos): [13, 14, 15, 16, 17, 18]
   ✅ CORREÇÃO: 6 nós inativos fixados em T=20.0°C
   🚨 CORREÇÃO FINAL: 2 nós órfãos detectados: [8, 12]
   ✅ CORREÇÃO FINAL: 2 nós órfãos adicionais fixados
   ✅ CORREÇÃO: 5 contornos ativos aplicados
   🎯 TOTAL FINAL: 13 condições de contorno (cobertura completa)
   ⚠️  GMRES falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver mais robusto (CG+ILU)...
   ⚠️  CG falhou: error code 73
[0] KSPSolve() at ./src/ksp/ksp/interface/itfunc.c:1070
[0] KSPSolve_Private() at ./src/ksp/ksp/interface/itfunc.c:824
[0] KSPSetUp() at ./src/ksp/ksp/interface/itfunc.c:406
[0] PCSetUp() at ./src/ksp/pc/interface/precon.c:994
[0] PCSetUp_ILU() at ./src/ksp/pc/impls/factor/ilu/ilu.c:135
[0] MatILUFactorSymbolic() at ./src/mat/interface/matrix.c:6897
[0] MatILUFactorSymbolic_SeqAIJ() at ./src/mat/impls/aij/seq/aijfact.c:1673
[0] Object is in wrong state
[0] Matrix is missing diagonal entry 13
   🔄 Tentando solver direto (LU)...
   ❌ Erro no solver de energia: Solução LU divergente
   🔄 Mantendo solução anterior (Tp não será atualizado)

✅ SIMULAÇÃO DO PRIMEIRO BLOCO CONCLUÍDA
   📊 Total: 22 passos
   🕒 Tempo final: 48.0h

🎯 SIMULAÇÃO PRIMEIRO BLOCO CONCLUÍDA
📝 Melhorias implementadas:
   ✅ Lógica assertiva usando JSON
   ✅ Formulação variacional apenas para domínios ativos
   ✅ Condições de contorno apenas para contornos ativos
   ✅ Execução apenas no primeiro bloco de tempo
   ✅ Função encontrar_bloco_para_tempo_atual implementada
