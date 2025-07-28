# REQUISITOS PARA REPLICAR O PROJETO FENICSx - BARRAGEM

## 📋 VISÃO GERAL

Este documento lista **todas as dependências** necessárias para replicar o projeto de simulação de barragens com FEniCSx em uma nova máquina.

---

## 🖥️ SISTEMA OPERACIONAL

### **Recomendado:**
- **Ubuntu 24.04 LTS** (ou versões mais recentes)
- **Linux** (qualquer distribuição baseada em Debian)

### **Alternativas:**
- **macOS** (com Homebrew)
- **Windows** (com WSL2 + Ubuntu)

---

## 🐍 PYTHON

### **Versão:**
- **Python 3.12+** (recomendado)
- **Python 3.11+** (mínimo)

### **Verificação:**
```bash
python --version
```

---

## 🔧 COMPILADOR C/C++

### **Requisitos:**
- **GCC 13.3.0+** (ou versão compatível)
- **G++** (compilador C++)
- **Make** (sistema de build)

### **Instalação (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install build-essential gcc g++ make
```

### **Verificação:**
```bash
gcc --version
g++ --version
make --version
```

---

## 📦 MPI (MESSAGE PASSING INTERFACE)

### **Requisitos:**
- **Open MPI 4.1.6+** (ou MPICH)
- **mpirun** (executor MPI)

### **Instalação (Ubuntu/Debian):**
```bash
sudo apt install openmpi-bin libopenmpi-dev
```

### **Verificação:**
```bash
mpirun --version
```

---

## 🐍 AMBIENTE VIRTUAL PYTHON

### **Criação:**
```bash
# Criar ambiente virtual
python -m venv fenics_env

# Ativar ambiente
source fenics_env/bin/activate  # Linux/macOS
# ou
fenics_env\Scripts\activate     # Windows
```

---

## 📚 DEPENDÊNCIAS PYTHON PRINCIPAIS

### **FEniCSx (Core):**
```bash
pip install fenics-dolfinx==0.9.0
pip install fenics-basix==0.9.0
pip install fenics-ffcx==0.9.0
pip install fenics-ufl==2024.2.0
```

### **MPI para Python:**
```bash
pip install mpi4py==3.1.5
```

### **PETSc (Solver):**
```bash
pip install petsc4py==3.19.6
pip install slepc4py==3.19.2
```

### **Processamento de Dados:**
```bash
pip install numpy==1.26.4
pip install scipy==1.11.4
pip install pandas==2.1.4
```

### **Visualização e I/O:**
```bash
pip install matplotlib==3.6.3
pip install h5py==3.10.0
pip install meshio==5.3.5
pip install pyvista==0.44.1
```

### **Configuração e Utilitários:**
```bash
pip install PyYAML==6.0.1
pip install rich==13.7.1
pip install tqdm
```

---

## 🛠️ FERRAMENTAS AUXILIARES

### **Gmsh (Gerador de Malhas):**
```bash
# Ubuntu/Debian
sudo apt install gmsh

# Verificação
gmsh --version
```

### **Git (Controle de Versão):**
```bash
sudo apt install git
```

### **Editor de Código (Opcional):**
```bash
# VS Code
sudo snap install code --classic

# ou Vim
sudo apt install vim
```

---

## 📋 SCRIPT DE INSTALAÇÃO COMPLETA

### **Criar arquivo `install_dependencies.sh`:**

```bash
#!/bin/bash

echo "🚀 INSTALANDO DEPENDÊNCIAS DO PROJETO FENICSx - BARRAGEM"
echo "========================================================"

# Atualizar sistema
echo "📦 Atualizando sistema..."
sudo apt update && sudo apt upgrade -y

# Instalar compiladores
echo "🔧 Instalando compiladores..."
sudo apt install -y build-essential gcc g++ make

# Instalar MPI
echo "📡 Instalando MPI..."
sudo apt install -y openmpi-bin libopenmpi-dev

# Instalar Gmsh
echo "📐 Instalando Gmsh..."
sudo apt install -y gmsh

# Instalar Git
echo "📝 Instalando Git..."
sudo apt install -y git

# Criar ambiente virtual
echo "🐍 Criando ambiente virtual Python..."
python3 -m venv fenics_env

# Ativar ambiente
echo "🔌 Ativando ambiente virtual..."
source fenics_env/bin/activate

# Atualizar pip
echo "📦 Atualizando pip..."
pip install --upgrade pip

# Instalar dependências Python
echo "📚 Instalando dependências Python..."

# FEniCSx Core
pip install fenics-dolfinx==0.9.0
pip install fenics-basix==0.9.0
pip install fenics-ffcx==0.9.0
pip install fenics-ufl==2024.2.0

# MPI e Solvers
pip install mpi4py==3.1.5
pip install petsc4py==3.19.6
pip install slepc4py==3.19.2

# Processamento de Dados
pip install numpy==1.26.4
pip install scipy==1.11.4
pip install pandas==2.1.4

# Visualização e I/O
pip install matplotlib==3.6.3
pip install h5py==3.10.0
pip install meshio==5.3.5
pip install pyvista==0.44.1

# Configuração e Utilitários
pip install PyYAML==6.0.1
pip install rich==13.7.1
pip install tqdm

echo "✅ INSTALAÇÃO CONCLUÍDA!"
echo "🎯 Para ativar o ambiente: source fenics_env/bin/activate"
```

### **Executar script:**
```bash
chmod +x install_dependencies.sh
./install_dependencies.sh
```

---

## 🔍 VERIFICAÇÃO DA INSTALAÇÃO

### **Script de verificação `verify_installation.py`:**

```python
#!/usr/bin/env python3
"""
Script para verificar se todas as dependências estão instaladas corretamente
"""

import sys
import importlib

def check_package(package_name, version=None):
    """Verifica se um pacote está instalado"""
    try:
        module = importlib.import_module(package_name)
        if version:
            if hasattr(module, '__version__'):
                print(f"✅ {package_name} {module.__version__}")
            else:
                print(f"✅ {package_name} (versão não disponível)")
        else:
            print(f"✅ {package_name}")
        return True
    except ImportError:
        print(f"❌ {package_name} NÃO INSTALADO")
        return False

def main():
    print("🔍 VERIFICANDO DEPENDÊNCIAS DO PROJETO FENICSx")
    print("=" * 50)
    
    # Verificar versão do Python
    print(f"🐍 Python {sys.version}")
    
    # Dependências principais
    packages = [
        ("dolfinx", "fenics-dolfinx"),
        ("basix", "fenics-basix"),
        ("ffcx", "fenics-ffcx"),
        ("ufl", "fenics-ufl"),
        ("mpi4py", "mpi4py"),
        ("petsc4py", "petsc4py"),
        ("slepc4py", "slepc4py"),
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("pandas", "pandas"),
        ("matplotlib", "matplotlib"),
        ("h5py", "h5py"),
        ("meshio", "meshio"),
        ("pyvista", "pyvista"),
        ("yaml", "PyYAML"),
        ("rich", "rich"),
        ("tqdm", "tqdm"),
    ]
    
    all_installed = True
    for module_name, package_name in packages:
        if not check_package(module_name):
            all_installed = False
    
    print("\n" + "=" * 50)
    if all_installed:
        print("🎉 TODAS AS DEPENDÊNCIAS ESTÃO INSTALADAS!")
        print("✅ O projeto está pronto para uso.")
    else:
        print("⚠️  ALGUMAS DEPENDÊNCIAS ESTÃO FALTANDO!")
        print("❌ Execute o script de instalação novamente.")

if __name__ == "__main__":
    main()
```

---

## 🚀 PRIMEIROS PASSOS

### **1. Clonar o repositório:**
```bash
git clone <URL_DO_REPOSITORIO>
cd FENICS
```

### **2. Ativar ambiente virtual:**
```bash
source fenics_env/bin/activate
```

### **3. Verificar instalação:**
```bash
python verify_installation.py
```

### **4. Testar execução:**
```bash
cd CalcT
python barragem-Gemini-R2.py barragem1
```

---

## 📊 VERSÕES TESTADAS

### **Sistema:**
- **OS:** Ubuntu 24.04 LTS
- **Python:** 3.12.3
- **GCC:** 13.3.0
- **Open MPI:** 4.1.6

### **Dependências Principais:**
- **FEniCSx:** 0.9.0
- **NumPy:** 1.26.4
- **SciPy:** 1.11.4
- **MPI4Py:** 3.1.5
- **PETSc4Py:** 3.19.6

---

## ⚠️ PROBLEMAS COMUNS

### **1. Erro de compilação:**
```bash
# Instalar dependências de desenvolvimento
sudo apt install python3-dev
```

### **2. Erro de MPI:**
```bash
# Verificar se MPI está funcionando
mpirun -np 2 python -c "from mpi4py import MPI; print('MPI OK')"
```

### **3. Erro de FEniCSx:**
```bash
# Reinstalar FEniCSx
pip uninstall fenics-dolfinx fenics-basix fenics-ffcx fenics-ufl
pip install fenics-dolfinx==0.9.0 fenics-basix==0.9.0 fenics-ffcx==0.9.0 fenics-ufl==2024.2.0
```

---

## 📞 SUPORTE

Se encontrar problemas durante a instalação:

1. **Verifique** se todas as versões estão corretas
2. **Consulte** a documentação oficial do FEniCSx
3. **Reporte** problemas específicos com logs detalhados

---

*Documento gerado em: 2025-07-28* 