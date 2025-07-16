#!/usr/bin/python3
"""
Teste automatizado para ETAPA 5 SIMPLIFICADA - Simulação Stagewise com Submeshes
"""

import subprocess
import sys
import os

def test_etapa5_simplificado():
    """Testa se a ETAPA 5 SIMPLIFICADA executa com sucesso"""
    print("=== TESTE ETAPA 5 SIMPLIFICADA - SIMULAÇÃO STAGEWISE ===")
    
    # Verificar se arquivos necessários existem
    required_files = [
        "Uploads/barragem2.xdmf",
        "Uploads/analise_stagewise_barragem2_xdmf.json",
        "etapa5_stagewise_submesh_simplificado.py"
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"✗ Arquivo necessário não encontrado: {file_path}")
            return False
        else:
            print(f"✓ Arquivo encontrado: {file_path}")
    
    try:
        # Executar ETAPA 5 SIMPLIFICADA
        print("\nExecutando ETAPA 5 SIMPLIFICADA...")
        result = subprocess.run([
            "python3", "etapa5_stagewise_submesh_simplificado.py"
        ], capture_output=True, text=True, timeout=900)  # 15 min timeout
        
        # Mostrar saída
        print("=== SAÍDA DO PROGRAMA ===")
        print(result.stdout)
        
        if result.stderr:
            print("=== ERROS ===")
            print(result.stderr)
        
        # Verificar código de saída
        if result.returncode == 0:
            print("✓ ETAPA 5 SIMPLIFICADA executou com sucesso")
            
            # Verificar se contém indicadores de sucesso
            if "RESULTADO: SUCESSO" in result.stdout:
                print("✓ Resultado final: SUCESSO")
            else:
                print("✗ Resultado final não encontrado")
                return False
                
            # Verificar se simulação stagewise funcionou
            if "SIMULANDO BLOCO 1" in result.stdout and \
               "SIMULANDO BLOCO 2" in result.stdout and \
               "SIMULANDO BLOCO 3" in result.stdout:
                print("✓ Simulação stagewise executada (3 blocos)")
            else:
                print("✗ Simulação stagewise não detectada")
                return False
                
            # Verificar se submeshes foram criadas
            if "Submesh criada:" in result.stdout:
                print("✓ Submeshes criadas com sucesso")
            else:
                print("✗ Criação de submeshes não detectada")
                return False
                
            # Verificar se temperatura evoluiu adequadamente
            if "Temperatura dentro de limites razoáveis" in result.stdout:
                print("✓ Temperatura permaneceu em limites razoáveis")
            else:
                print("✗ Possível divergência de temperatura")
                return False
                
            # Verificar se hidratação progrediu
            if "Hidratação progrediu adequadamente" in result.stdout:
                print("✓ Hidratação progrediu adequadamente")
            else:
                print("✗ Hidratação não progrediu adequadamente")
                return False
                
            # Verificar eficiência stagewise
            if "SIMULAÇÃO STAGEWISE" in result.stdout:
                print("✓ Simulação stagewise detectada")
            else:
                print("✗ Simulação stagewise não detectada")
                return False
                
            # Verificar se relatório final está completo
            checks = [
                "Simulação stagewise concluída",
                "Temperatura inicial:",
                "Temperatura final:",
                "Pico de temperatura:",
                "Hidratação final:"
            ]
            
            for check in checks:
                if check in result.stdout:
                    print(f"✓ {check} - OK")
                else:
                    print(f"✗ {check} - NÃO ENCONTRADO")
                    return False
                
            return True
        else:
            print(f"✗ ETAPA 5 SIMPLIFICADA falhou com código {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 5 SIMPLIFICADA excedeu tempo limite (15 min)")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 5 SIMPLIFICADA: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa5_simplificado()
    if success:
        print("\nTESTE ETAPA 5 SIMPLIFICADA: PASSOU")
        print("🎉 Simulação stagewise com submeshes FUNCIONANDO!")
        print("💡 Implementação básica de construção por etapas concluída")
        sys.exit(0)
    else:
        print("\nTESTE ETAPA 5 SIMPLIFICADA: FALHOU")
        sys.exit(1) 