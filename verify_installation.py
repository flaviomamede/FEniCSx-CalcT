#!/usr/bin/env python3
"""
Script para verificar se todas as dependências estão instaladas corretamente
"""

import sys
import importlib
import subprocess

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

def check_system_tool(tool_name, version_flag="--version"):
    """Verifica se uma ferramenta do sistema está instalada"""
    try:
        result = subprocess.run([tool_name, version_flag], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version_line = result.stdout.strip().split('\n')[0]
            print(f"✅ {tool_name}: {version_line}")
            return True
        else:
            print(f"❌ {tool_name} NÃO FUNCIONANDO")
            return False
    except FileNotFoundError:
        print(f"❌ {tool_name} NÃO INSTALADO")
        return False
    except subprocess.TimeoutExpired:
        print(f"⚠️  {tool_name} TIMEOUT")
        return False

def main():
    print("🔍 VERIFICANDO DEPENDÊNCIAS DO PROJETO FENICSx")
    print("=" * 60)
    
    # Verificar versão do Python
    print(f"🐍 Python {sys.version}")
    print()
    
    # Verificar ferramentas do sistema
    print("🛠️  FERRAMENTAS DO SISTEMA:")
    print("-" * 30)
    system_tools = [
        ("gcc", "--version"),
        ("g++", "--version"),
        ("make", "--version"),
        ("mpirun", "--version"),
        ("gmsh", "--version"),
        ("git", "--version"),
    ]
    
    system_ok = True
    for tool, flag in system_tools:
        if not check_system_tool(tool, flag):
            system_ok = False
    
    print()
    
    # Verificar dependências Python
    print("📚 DEPENDÊNCIAS PYTHON:")
    print("-" * 30)
    
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
    
    python_ok = True
    for module_name, package_name in packages:
        if not check_package(module_name):
            python_ok = False
    
    print()
    print("=" * 60)
    
    # Resumo final
    if system_ok and python_ok:
        print("🎉 TODAS AS DEPENDÊNCIAS ESTÃO INSTALADAS!")
        print("✅ O projeto está pronto para uso.")
        print()
        print("🚀 PRÓXIMOS PASSOS:")
        print("1. Ativar ambiente: source fenics_env/bin/activate")
        print("2. Testar projeto: cd CalcT && python barragem-Gemini-R2.py barragem1")
    else:
        print("⚠️  ALGUMAS DEPENDÊNCIAS ESTÃO FALTANDO!")
        print("❌ Execute o script de instalação novamente:")
        print("   chmod +x install_dependencies.sh && ./install_dependencies.sh")
    
    print()
    print("📊 RESUMO:")
    print(f"   Sistema: {'✅' if system_ok else '❌'}")
    print(f"   Python:  {'✅' if python_ok else '❌'}")
    print(f"   Geral:   {'✅' if (system_ok and python_ok) else '❌'}")

if __name__ == "__main__":
    main() 