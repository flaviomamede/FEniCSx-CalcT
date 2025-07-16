#!/usr/bin/python3
"""
Teste automatizado para ETAPA 3 - Sistema Completo (Versão Local)
"""

import subprocess
import sys
import os

def test_etapa3_local():
    """Testa se a ETAPA 3 executa com sucesso no ambiente local"""
    print("=== TESTE ETAPA 3 LOCAL ===")
    
    # Verificar se arquivos necessários existem (caminhos relativos ao diretório atual)
    required_files = [
        "Uploads/barragem1.xdmf",
        "Uploads/analise_stagewise_barragem1_xdmf.json",
        "etapa3_sistema_completo_local.py"
    ]
    
    for file_path in required_files:
        if not os.path.exists(file_path):
            print(f"✗ Arquivo necessário não encontrado: {file_path}")
            return False
        else:
            print(f"✓ Arquivo encontrado: {file_path}")
    
    try:
        # Executar ETAPA 3
        print("\nExecutando ETAPA 3...")
        result = subprocess.run([
            "python3", "etapa3_sistema_completo_local.py"
        ], capture_output=True, text=True, timeout=300)
        
        # Mostrar saída
        print("=== SAÍDA DO PROGRAMA ===")
        print(result.stdout)
        
        if result.stderr:
            print("=== ERROS ===")
            print(result.stderr)
        
        # Verificar código de saída
        if result.returncode == 0:
            print("✓ ETAPA 3 executou com sucesso")
            
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
            
            # Verificar se arquivos de saída foram criados
            output_files = [
                "etapa3_temperatura.xdmf",
                "etapa3_teq.xdmf",
                "etapa3_Q.xdmf"
            ]
            
            for output_file in output_files:
                if os.path.exists(output_file):
                    print(f"✓ Arquivo de saída criado: {output_file}")
                else:
                    print(f"✗ Arquivo de saída não encontrado: {output_file}")
                    return False
                
            return True
        else:
            print(f"✗ ETAPA 3 falhou com código {result.returncode}")
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 3 excedeu tempo limite")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 3: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa3_local()
    if success:
        print("\nTESTE ETAPA 3: PASSOU")
        sys.exit(0)
    else:
        print("\nTESTE ETAPA 3: FALHOU")
        sys.exit(1)