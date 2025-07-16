
#!/usr/bin/env python3
"""
Script de teste para validar os resultados da ETAPA 1
"""

import numpy as np
import subprocess
import sys
import os

def test_etapa1():
    """Testa se a ETAPA 1 executou corretamente"""
    
    print("=== TESTE DA ETAPA 1 ===")
    
    # 1. Executar o script principal
    print("1. Executando etapa1_conducao_pura.py...")
    
    try:
        result = subprocess.run([
            sys.executable, "/home/ubuntu/etapa1_conducao_pura.py"
        ], capture_output=True, text=True, timeout=300)
        
        print(f"   Código de saída: {result.returncode}")
        
        if result.returncode != 0:
            print("   ✗ Script falhou na execução")
            print("   STDOUT:", result.stdout[-500:])  # Últimas 500 chars
            print("   STDERR:", result.stderr[-500:])
            return False
        else:
            print("   ✓ Script executou sem erros")
            
    except subprocess.TimeoutExpired:
        print("   ✗ Script excedeu timeout de 300s")
        return False
    except Exception as e:
        print(f"   ✗ Erro na execução: {e}")
        return False
    
    # 2. Verificar se não houve erro "Matrix is missing diagonal entry"
    print("2. Verificando erros de matriz...")
    
    output_text = result.stdout + result.stderr
    if "Matrix is missing diagonal entry" in output_text:
        print("   ✗ Erro 'Matrix is missing diagonal entry' detectado!")
        return False
    else:
        print("   ✓ Nenhum erro de matriz diagonal detectado")
    
    # 3. Verificar se arquivo de saída foi criado
    print("3. Verificando arquivo de saída...")
    
    output_file = "/home/ubuntu/etapa1_resultado.xdmf"
    if os.path.exists(output_file):
        print(f"   ✓ Arquivo {output_file} criado")
    else:
        print(f"   ✗ Arquivo {output_file} não encontrado")
        return False
    
    # 4. Verificar se solução está próxima de 20°C
    print("4. Verificando valores de temperatura...")
    
    # Procurar por valores de temperatura no output
    lines = result.stdout.split('\n')
    temp_values = []
    
    for line in lines:
        if "T_mean =" in line:
            try:
                temp_str = line.split("=")[1].split("°C")[0].strip()
                temp_val = float(temp_str)
                temp_values.append(temp_val)
            except:
                pass
    
    if temp_values:
        final_temp = temp_values[-1]  # Última temperatura média
        print(f"   Temperatura final média: {final_temp:.6f} °C")
        
        # Verificar se está próxima de 20°C (tolerância de 5°C)
        if abs(final_temp - 20.0) < 5.0:
            print("   ✓ Temperatura dentro da faixa esperada (15-25°C)")
        else:
            print("   ✗ Temperatura fora da faixa esperada")
            return False
    else:
        print("   ⚠ Não foi possível extrair valores de temperatura")
    
    # 5. Verificar se houve convergência
    print("5. Verificando convergência...")
    
    if "ETAPA 1 CONCLUÍDA COM SUCESSO" in result.stdout:
        print("   ✓ Simulação convergiu com sucesso")
    else:
        print("   ✗ Simulação não convergiu")
        return False
    
    print("\n=== TODOS OS TESTES PASSARAM ===")
    return True

if __name__ == "__main__":
    success = test_etapa1()
    if success:
        print("RESULTADO FINAL: ✓ ETAPA 1 VALIDADA")
        exit(0)
    else:
        print("RESULTADO FINAL: ✗ ETAPA 1 FALHOU")
        exit(1)
