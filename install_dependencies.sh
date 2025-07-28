#!/bin/bash

echo "🚀 INSTALANDO DEPENDÊNCIAS DO PROJETO FENICSx - BARRAGEM"
echo "========================================================"

# Atualizar sistema
echo "📦 Atualizando sistema..."
sudo apt update && sudo apt upgrade -y

# Instalar compiladores
echo "🔧 Instalando compiladores..."
sudo apt install -y build-essential gcc g++ make python3-dev

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
echo "  - Instalando FEniCSx Core..."
pip install fenics-dolfinx==0.9.0
pip install fenics-basix==0.9.0
pip install fenics-ffcx==0.9.0
pip install fenics-ufl==2024.2.0

# MPI e Solvers
echo "  - Instalando MPI e Solvers..."
pip install mpi4py==3.1.5
pip install petsc4py==3.19.6
pip install slepc4py==3.19.2

# Processamento de Dados
echo "  - Instalando bibliotecas de processamento..."
pip install numpy==1.26.4
pip install scipy==1.11.4
pip install pandas==2.1.4

# Visualização e I/O
echo "  - Instalando bibliotecas de visualização..."
pip install matplotlib==3.6.3
pip install h5py==3.10.0
pip install meshio==5.3.5
pip install pyvista==0.44.1

# Configuração e Utilitários
echo "  - Instalando utilitários..."
pip install PyYAML==6.0.1
pip install rich==13.7.1
pip install tqdm

echo ""
echo "✅ INSTALAÇÃO CONCLUÍDA!"
echo "🎯 Para ativar o ambiente: source fenics_env/bin/activate"
echo ""
echo "🔍 Para verificar a instalação: python verify_installation.py"
echo "🚀 Para testar o projeto: cd CalcT && python barragem-Gemini-R2.py barragem1" 