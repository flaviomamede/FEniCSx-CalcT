#!/usr/bin/python3

"""
Teste da Etapa 5 - Simulau00e7u00e3o Stagewise com Interpolau00e7u00e3o Correta
"""

import os
import sys

def verificar_arquivos():
    """Verifica se os arquivos necessu00e1rios existem"""
    arquivos_necessarios = [
        "Uploads/barragem2.xdmf",
        "Uploads/analise_stagewise_barragem2_xdmf.json",
        "etapa5_stagewise_submesh_copy.py"
    ]
    
    for arquivo in arquivos_necessarios:
        if os.path.exists(arquivo):
            print(f"u2713 Arquivo encontrado: {arquivo}")
        else:
            print(f"u2717 Arquivo nu00e3o encontrado: {arquivo}")
            return False
    
    return True

def executar_teste():
    """Executa o teste da Etapa 5"""
    print("\nExecutando ETAPA 5 com interpolau00e7u00e3o correta...")
    
    # Executar o script
    codigo_saida = os.system("python3 etapa5_stagewise_submesh_copy.py")
    
    if codigo_saida == 0:
        print(f"u2713 ETAPA 5 executou com sucesso")
        return True
    else:
        print(f"u2717 ETAPA 5 falhou com cu00f3digo {codigo_saida}")
        return False

def main():
    print("=== TESTE ETAPA 5 - INTERPOLAu00c7u00c3O CORRETA ====")
    
    # Verificar arquivos
    if not verificar_arquivos():
        print("\nTESTE ETAPA 5: FALHOU (arquivos nu00e3o encontrados)")
        return False
    
    # Executar teste
    if not executar_teste():
        print("\nTESTE ETAPA 5: FALHOU (execuu00e7u00e3o)")
        return False
    
    print("\nTESTE ETAPA 5: PASSOU")
    return True

if __name__ == "__main__":
    sucesso = main()
    sys.exit(0 if sucesso else 1)