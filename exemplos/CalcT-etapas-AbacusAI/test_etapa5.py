#!/usr/bin/python3
"""
Teste automatizado para ETAPA 5 - Simulação Stagewise com Interpolação Correta
"""

import subprocess
import sys
import os

def test_etapa5_interpolacao():
    """Testa se a ETAPA 5 executa com interpolação correta entre blocos"""
    print("=== TESTE ETAPA 5 - INTERPOLAÇÃO CORRETA ===")
    
    # Verificar se arquivos necessários existem
    required_files = [
        "Uploads/barragem2.xdmf",
        "Uploads/barragem2.yaml",
        "Uploads/analise_stagewise_barragem2_xdmf.json",
        "etapa5_stagewise_submesh.py"
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"✗ Arquivo necessário não encontrado: {file_path}")
            return False
        else:
            print(f"✓ Arquivo encontrado: {file_path}")
    
    try:
        # Executar ETAPA 5
        print("\nExecutando ETAPA 5 com interpolação correta...")
        result = subprocess.run([
            "python3", "etapa5_stagewise_submesh.py"
        ], capture_output=True, text=True, timeout=600)
        
        # Mostrar saída
        print("=== SAÍDA DO PROGRAMA ===")
        print(result.stdout)
        
        if result.stderr:
            print("=== ERROS ===")
            print(result.stderr)
        
        # Verificar código de saída
        if result.returncode == 0:
            print("✓ ETAPA 5 executou com sucesso")
            
            # Verificar se contém indicadores de sucesso
            if "RESULTADO: SUCESSO" in result.stdout:
                print("✓ Resultado final: SUCESSO")
            else:
                print("✗ Resultado final não encontrado")
                return False
            
            # Verificar se interpolação foi executada
            if "PASSO 3.1: Interpolando estado do bloco anterior" in result.stdout:
                print("✓ Interpolação entre blocos executada")
            else:
                print("✗ Interpolação entre blocos não detectada")
                return False
            
            # Verificar se aplicação de C.I. nos novos elementos foi executada
            if "PASSO 3.2: Aplicando C.I. nos elementos novos" in result.stdout:
                print("✓ Condições iniciais aplicadas nos novos elementos")
            else:
                print("✗ Condições iniciais nos novos elementos não detectadas")
                return False
            
            # Verificar se DOFs novos foram identificados
            if "DOFs novos (T):" in result.stdout:
                print("✓ DOFs novos identificados corretamente")
            else:
                print("✗ DOFs novos não identificados")
                return False
            
            # Verificar se sistema acoplado funcionou
            if "SIMULAÇÃO STAGEWISE COMPLETA COM INTERPOLAÇÃO CORRETA" in result.stdout:
                print("✓ Simulação stagewise com interpolação completa")
            else:
                print("✗ Simulação stagewise não concluída")
                return False
            
            # Verificar se não houve segfault
            if "Segmentation fault" not in result.stderr and "segfault" not in result.stderr.lower():
                print("✓ Nenhum segmentation fault detectado")
            else:
                print("✗ Segmentation fault detectado")
                return False
            
            # Verificar se todos os 3 blocos foram processados
            blocos_encontrados = result.stdout.count("BLOCO") - result.stdout.count("BLOCO CONSTRUTIVO")
            if blocos_encontrados >= 3:
                print(f"✓ Todos os blocos processados: {blocos_encontrados}")
            else:
                print(f"✗ Apenas {blocos_encontrados} blocos processados")
                return False
            
            # Verificar se submeshes foram criadas
            if "Submesh criada:" in result.stdout:
                print("✓ Submeshes criadas corretamente")
            else:
                print("✗ Submeshes não criadas")
                return False
            
            # Verificar se elementos cumulativos foram processados
            if "Elementos ativos totais:" in result.stdout:
                print("✓ Elementos cumulativos processados")
            else:
                print("✗ Elementos cumulativos não processados")
                return False
            
            return True
        else:
            print(f"✗ ETAPA 5 falhou com código {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 5 excedeu tempo limite")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 5: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa5_interpolacao()
    if success:
        print("\n🎉 TESTE ETAPA 5: PASSOU - INTERPOLAÇÃO CORRETA IMPLEMENTADA!")
        sys.exit(0)
    else:
        print("\n❌ TESTE ETAPA 5: FALHOU - PROBLEMAS NA INTERPOLAÇÃO")
        sys.exit(1) 