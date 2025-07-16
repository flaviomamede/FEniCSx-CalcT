
#!/usr/bin/python3
"""
Teste automatizado para ETAPA 3 - Sistema Completo
"""

import subprocess
import sys

def test_etapa3():
    """Testa se a ETAPA 3 executa com sucesso"""
    print("=== TESTE ETAPA 3 ===")
    
    try:
        # Executar ETAPA 3
        result = subprocess.run(["/usr/bin/python3", "etapa3_sistema_completo.py"], 
                              capture_output=True, text=True, timeout=180)
        
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
                
            return True
        else:
            print(f"✗ ETAPA 3 falhou com código {result.returncode}")
            print("STDERR:", result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 3 excedeu tempo limite")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 3: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa3()
    if success:
        print("TESTE ETAPA 3: PASSOU")
        sys.exit(0)
    else:
        print("TESTE ETAPA 3: FALHOU")
        sys.exit(1)
