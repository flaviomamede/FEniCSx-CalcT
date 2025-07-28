# Guia de Uso - Salvamento de Submeshes

## Visão Geral

Agora o sistema salva automaticamente todas as submeshes criadas para cada camada na pasta do caso, facilitando análise posterior e validação dos resultados.

## Estrutura de Arquivos

### Pasta de Saída
```
caso/
├── submeshes/
│   ├── caso_camada_1.xdmf
│   ├── caso_camada_1.h5
│   ├── caso_camada_2.xdmf
│   ├── caso_camada_2.h5
│   └── results/
│       ├── caso_camada_1_step0000.xdmf
│       ├── caso_camada_1_step0000.h5
│       ├── caso_camada_2_step0000.xdmf
│       └── caso_camada_2_step0000.h5
├── resultados/
└── caso_resultados.xdmf
```

## Como Usar

### 1. Ativar Criação de Submeshes

Para usar a funcionalidade de submeshes, você precisa primeiro criar as submeshes:

```python
# No seu código principal
simulator = BarragemFEniCSx(config_file="seu_caso.yaml")

# Criar submeshes para cada camada
simulator._create_layer_submeshes()

# Executar simulação
simulator.run()
```

### 2. Arquivos Gerados

#### Submeshes Estruturais
- **caso_camada_N.xdmf**: Malha da camada N
- **caso_camada_N.h5**: Dados da malha (HDF5)

#### Resultados por Camada
- **caso_camada_N_stepMMMM.xdmf**: Resultados da camada N no passo MMMM
- **caso_camada_N_stepMMMM.h5**: Dados dos resultados (HDF5)

### 3. Configuração YAML

Não é necessário nenhuma configuração adicional no YAML. O sistema salva automaticamente quando as submeshes são criadas.

## Exemplo de Uso Completo

```python
#!/usr/bin/env python3

from pathlib import Path
import sys
import os

# Adicionar o diretório atual ao path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from barragem-Gemini-R2 import BarragemFEniCSx

def main():
    # Caminho para o arquivo de configuração
    config_file = "barragem1/barragem1.yaml"
    
    # Criar simulador
    simulator = BarragemFEniCSx(config_file=config_file)
    
    # Criar submeshes (opcional - se quiser usar submeshes)
    simulator._create_layer_submeshes()
    
    # Executar simulação
    simulator.run()
    
    print("\n" + "="*60)
    print("SIMULAÇÃO CONCLUÍDA!")
    print("="*60)
    
    # Verificar arquivos salvos
    case_dir = Path(config_file).parent
    submesh_dir = case_dir / "submeshes"
    
    if submesh_dir.exists():
        print(f"\nSubmeshes salvas em: {submesh_dir}")
        for file in submesh_dir.glob("*.xdmf"):
            print(f"   - {file.name}")
        
        results_dir = submesh_dir / "results"
        if results_dir.exists():
            print(f"\nResultados por camada em: {results_dir}")
            for file in results_dir.glob("*.xdmf"):
                print(f"   - {file.name}")

if __name__ == "__main__":
    main()
```

## Visualização das Submeshes

### Paraview
```bash
# Abrir submesh estrutural
paraview caso/submeshes/caso_camada_1.xdmf

# Abrir resultados por camada
paraview caso/submeshes/results/caso_camada_1_step0000.xdmf
```

### Python com PyVista
```python
import pyvista as pv
import dolfinx.plot
from dolfinx import io

# Carregar submesh
case_dir = Path("barragem1")
submesh_file = case_dir / "submeshes" / "barragem1_camada_1.xdmf"

# Ler malha
mesh = pv.read(str(submesh_file))

# Visualizar
plotter = pv.Plotter()
plotter.add_mesh(mesh, show_edges=True)
plotter.show()
```

## Scripts Úteis

### 1. Verificar Submeshes Criadas
```python
#!/usr/bin/env python3
# check_submeshes.py

import os
from pathlib import Path

def check_submeshes(caso_path):
    """Verifica quais submeshes foram criadas para um caso."""
    
    caso_path = Path(caso_path)
    submesh_dir = caso_path / "submeshes"
    
    if not submesh_dir.exists():
        print(f"Nenhuma submesh encontrada em {caso_path}")
        return
    
    print(f"\nSubmeshes encontradas em {caso_path}:")
    
    # Listar submeshes estruturais
    xdmf_files = list(submesh_dir.glob("*.xdmf"))
    if xdmf_files:
        print("\nEstruturais:")
        for file in sorted(xdmf_files):
            print(f"   ✅ {file.name}")
    
    # Listar resultados
    results_dir = submesh_dir / "results"
    if results_dir.exists():
        result_files = list(results_dir.glob("*.xdmf"))
        if result_files:
            print("\nResultados:")
            for file in sorted(result_files):
                print(f"   📊 {file.name}")
    
    # Estatísticas
    total_files = len(xdmf_files) + len(result_files) if results_dir.exists() else len(xdmf_files)
    print(f"\nTotal de arquivos: {total_files}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        check_submeshes(sys.argv[1])
    else:
        check_submeshes("barragem1")
```

### 2. Comparar Submeshes
```python
#!/usr/bin/env python3
# compare_submeshes.py

import numpy as np
from dolfinx import io
from pathlib import Path

def compare_submeshes(caso_path):
    """Compara informações entre submeshes."""
    
    caso_path = Path(caso_path)
    submesh_dir = caso_path / "submeshes"
    
    if not submesh_dir.exists():
        print("Nenhuma submesh encontrada")
        return
    
    print(f"\nComparação de Submeshes - {caso_path.name}")
    print("=" * 50)
    
    for xdmf_file in sorted(submesh_dir.glob("*.xdmf")):
        layer_name = xdmf_file.stem.replace(f"{caso_path.name}_", "")
        
        # Ler informações da submesh
        with io.XDMFFile(MPI.COMM_WORLD, str(xdmf_file), "r") as xdmf:
            mesh = xdmf.read_mesh(name=layer_name)
            
            # Estatísticas
            num_cells = mesh.topology.index_map(mesh.topology.dim).size_global
            num_vertices = mesh.topology.index_map(0).size_global
            
            print(f"\n{layer_name}:")
            print(f"   Células: {num_cells:,}")
            print(f"   Vértices: {num_vertices:,}")

if __name__ == "__main__":
    from mpi4py import MPI
    import sys
    
    if len(sys.argv) > 1:
        compare_submeshes(sys.argv[1])
    else:
        compare_submeshes("barragem1")
```

## Dicas de Uso

### 1. Economia de Espaço
As submeshes são salvas automaticamente, mas você pode controlar:

```python
# Desabilitar salvamento de resultados por camada (apenas estrutura)
simulator.save_submesh_results = False

# Salvar apenas em passos específicos
if step % 10 == 0:
    simulator._save_all_submeshes_and_results(step, current_time)
```

### 2. Integração com Scripts Existentes
```python
# Modificar seu script existente para incluir submeshes
# Basta adicionar uma linha após criar o simulador:

simulator = BarragemFEniCSx(config_file="seu_caso.yaml")
simulator._create_layer_submeshes()  # Adicionar esta linha
simulator.run()
```

### 3. Verificação Rápida
```bash
# Verificar se submeshes foram criadas
ls barragem1/submeshes/

# Verificar resultados por camada
ls barragem1/submeshes/results/
```

## Solução de Problemas

### Problema: Submeshes não aparecem
**Causa**: As submeshes não foram criadas
**Solução**: Adicionar `simulator._create_layer_submeshes()` antes de `simulator.run()`

### Problema: Arquivos muito grandes
**Causa**: Salvamento frequente de resultados
**Solução**: Ajustar frequência de salvamento no YAML

### Problema: Erro de permissão
**Causa**: Diretório não existe ou sem permissão
**Solução**: Verificar permissões e criar diretórios manualmente