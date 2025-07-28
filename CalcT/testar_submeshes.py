#!/usr/bin/env python3
"""
Script de teste para verificar salvamento de submeshes
"""

import sys
from pathlib import Path
from barragem_gemini_r2 import SimulacaoBarragemR2

def testar_submeshes(case_dir='barragem2'):
    """Testa criação e salvamento de submeshes"""
    print("=== TESTE DE SUBMESHES ===")

    case_path = Path(case_dir)
    config_file = case_path / f"{case_path.name}.yaml"
    json_file = case_path / f"{case_path.name}-msh.json"

    if not config_file.exists():
        print(f"Arquivo de configuração não encontrado: {config_file}")
        return False
    if not json_file.exists():
        print(f"Arquivo JSON da malha não encontrado: {json_file}")
        return False

    print(f"1. Configuração: {config_file}")
    print(f"2. JSON da malha: {json_file}")

    # Criar simulador
    simulator = SimulacaoBarragemR2(
        config_file=str(config_file),
        json_file=str(json_file)
    )
    
    print("2. Simulador criado")
    
    # Criar submeshes
    try:
        simulator._create_layer_submeshes()
        print(f"3. Submeshes criadas: {len(simulator.layer_submeshes)}")
        
        for layer_id, info in simulator.layer_submeshes.items():
            num_elements = len(info['active_elements'])
            print(f"   - Camada {layer_id+1}: {num_elements} elementos")
            
    except Exception as e:
        print(f"Erro ao criar submeshes: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Salvar submeshes
    try:
        simulator._save_submeshes(caso_path)
        print("4. Submeshes salvas!")
        
        # Verificar arquivos criados
        submesh_dir = caso_path / "submeshes"
        if submesh_dir.exists():
            files = list(submesh_dir.glob("*.xdmf"))
            print(f"5. Arquivos criados: {len(files)}")
            for f in sorted(files):
                print(f"   - {f.name}")
                
            # Verificar tamanho dos arquivos
            for f in files:
                h5_file = f.with_suffix('.h5')
                if h5_file.exists():
                    size = h5_file.stat().st_size / 1024  # KB
                    print(f"     {h5_file.name}: {size:.1f} KB")
        else:
            print("5. Diretório submeshes não foi criado")
            
    except Exception as e:
        print(f"Erro ao salvar submeshes: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    return True

if __name__ == "__main__":
    caso = sys.argv[1] if len(sys.argv) > 1 else 'barragem2'
    testar_submeshes(caso)