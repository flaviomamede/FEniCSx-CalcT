
#!/usr/bin/python3
"""
Teste automatizado para ETAPA 2 - Fonte Constante
"""

import subprocess
import sys

def test_etapa2():
    """Testa se a ETAPA 2 executa com sucesso"""
    print("=== TESTE ETAPA 2 ===")
    
    try:
        # Executar ETAPA 2
        result = subprocess.run(["/usr/bin/python3", "etapa2_fonte_constante.py"], 
                              capture_output=True, text=True, timeout=120)
        
        # Verificar código de saída
        if result.returncode == 0:
            print("✓ ETAPA 2 executou com sucesso")
            
            # Verificar se contém indicadores de sucesso
            if "RESULTADO: SUCESSO" in result.stdout:
                print("✓ Resultado final: SUCESSO")
            else:
                print("✗ Resultado final não encontrado")
                return False
                
            # Verificar se temperatura aumentou
            if "Temperatura aumentou devido à fonte de calor" in result.stdout:
                print("✓ Temperatura aumentou conforme esperado")
            else:
                print("✗ Aumento de temperatura não detectado")
                return False
                
            return True
        else:
            print(f"✗ ETAPA 2 falhou com código {result.returncode}")
            print("STDERR:", result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ ETAPA 2 excedeu tempo limite")
        return False
    except Exception as e:
        print(f"✗ Erro ao executar ETAPA 2: {e}")
        return False

if __name__ == "__main__":
    success = test_etapa2()
    if success:
        print("TESTE ETAPA 2: PASSOU")
        sys.exit(0)
    else:
        print("TESTE ETAPA 2: FALHOU")
        sys.exit(1)
