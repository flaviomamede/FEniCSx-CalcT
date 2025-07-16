#!/usr/bin/python3
"""
Teste automatizado para ETAPA 4 - Sistema com Submeshes
"""

import subprocess
import sys
import os

def test_etapa4():
    """Testa se a ETAPA 4 executa com sucesso"""
    print("=== TESTE ETAPA 4 - SUBMESHES ===")
    
    # Verificar se arquivos necessários existem (caminhos relativos ao diretório atual)
    required_files = [
        "Uploads/barragem2.xdmf",
        "Uploads/analise_stagewise_barragem2_xdmf.json",
        "etapa4_submesh.py"
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"✗ Arquivo necessário não encontrado: {file_path}")
            return False
        else:
            print(f"✓ Arquivo encontrado: {file_path}")
    
    try:
        # Executar ETAPA 4
        print("\nExecutando ETAPA 4...")
        result = subprocess.run([
            "python3", "etapa4_submesh.py"
        ], capture_output=True, text=True, timeout=300)
        
        # Mostrar saída
        print("=== SAÍDA DO PROGRAMA ===")
        print(result.stdout)
        
        if result.stderr:
            print("=== ERROS ===")
            print(result.stderr)
        
        # Verificar código de saída
        if result.returncode == 0:
            print("✓ ETAPA 4 executou com sucesso")
            
            # Verificar se contém indicadores de sucesso
            if "RESULTADO: SUCESSO" in result.stdout:
                print("✓ Resultado final: SUCESSO")
            else:
                print("✗ Resultado final não encontrado")
                return False
                
            # Verificar se sistema acoplado funcionou
            if "SIMULAÇÃO ACOPLADA" in result.stdout:
                print("✓ Simulação acoplada executada")
            else:
                print("✗ Simulação acoplada não detectada")
                return False
                
            # Verificar se não divergiu
            if "Temperatura dentro de limites razoáveis" in result.stdout:
                print("✓ Temperatura permaneceu em limites razoáveis")
            else:
                print("✗ Possível divergência de temperatura")
                return False
                
            # Verificar se tempo equivalente evoluiu
            if "teq_mean" in result.stdout:
                print("✓ Tempo equivalente calculado")
            else:
                print("✗ Tempo equivalente não encontrado")
                return False
                
            # Verificar se submesh foi criada corretamente
            if "Submesh criada: 54 células" in result.stdout:
                print("✓ Submesh criada corretamente")
            else:
                print("✗ Submesh não detectada")
                return False
                
            # Verificar se houve redução da malha
            if "Redução: 60.0% da malha global" in result.stdout:
                print("✓ Redução de malha confirmada")
            else:
                print("✗ Redução de malha não detectada")
                return False
            
            # Verificar se arquivos de saída foram criados (ETAPA 4)
            output_files = [
                "etapa4_temperatura.xdmf",
                "etapa4_teq.xdmf",
                "etapa4_Q.xdmf"
            ]
            
            for output_file in output_files:
                if os.path.exists(output_file):
                    print(f"✓ Arquivo de saída criado: {output_file}")
                else:
                    print(f"✗ Arquivo de saída não encontrado: {output_file}")
                    return False
                
            return True
        else:
            print(f"✗ ETAPA 4 falhou com código {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 4 excedeu tempo limite")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 4: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa4()
    if success:
        print("\nTESTE ETAPA 4: PASSOU")
        sys.exit(0)
    else:
        print("\nTESTE ETAPA 4: FALHOU")
        sys.exit(1) 