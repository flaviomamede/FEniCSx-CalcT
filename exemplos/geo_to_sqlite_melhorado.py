#!/usr/bin/env python3
"""
Conversor melhorado de arquivos .geo do GMSH para banco de dados SQLite
Versão segura e robusta com parser de expressões melhorado
"""

import sqlite3
import re
import ast
import operator
from typing import Dict, Any, List

class SafeExpressionEvaluator:
    """Avaliador seguro de expressões matemáticas simples"""
    
    # Operadores permitidos
    operators = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.USub: operator.neg,
    }
    
    def __init__(self, variables: Dict[str, float]):
        self.variables = variables
    
    def evaluate(self, expression: str) -> float:
        """Avalia uma expressão matemática de forma segura"""
        try:
            # Primeiro, substituir variáveis conhecidas
            expr = expression.strip()
            for var_name, value in self.variables.items():
                expr = expr.replace(var_name, str(value))
            
            # Se é apenas um número, retornar diretamente
            try:
                return float(expr)
            except ValueError:
                pass
            
            # Avaliar expressão usando AST (seguro)
            node = ast.parse(expr, mode='eval')
            return self._eval_node(node.body)
        
        except Exception as e:
            print(f"⚠️ Erro ao avaliar expressão '{expression}': {e}")
            return 0.0
    
    def _eval_node(self, node):
        """Avalia um nó AST recursivamente"""
        if isinstance(node, ast.Constant):  # Python 3.8+
            return node.value
        elif isinstance(node, ast.Num):  # Python < 3.8
            return node.n
        elif isinstance(node, ast.Name):
            if node.id in self.variables:
                return self.variables[node.id]
            else:
                raise ValueError(f"Variável desconhecida: {node.id}")
        elif isinstance(node, ast.BinOp):
            return self.operators[type(node.op)](
                self._eval_node(node.left),
                self._eval_node(node.right)
            )
        elif isinstance(node, ast.UnaryOp):
            return self.operators[type(node.op)](self._eval_node(node.operand))
        else:
            raise TypeError(f"Tipo de nó não suportado: {type(node)}")

def parse_geo_file(filepath: str) -> Dict[str, Any]:
    """
    Parser melhorado para arquivos .geo
    """
    print(f"📖 Analisando arquivo: {filepath}")
    
    data = {
        'variables': {},
        'points': [],
        'lines': [],
        'curve_loops': {},
        'plane_surfaces': [],
        'physical_surfaces': [],
        'physical_lines': [],
        'mesh_configs': {},
        'comments': []
    }

    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"❌ Erro ao ler arquivo: {e}")
        return data

    # 1. Parse Variables
    print("🔍 Extraindo variáveis...")
    variables_matches = re.findall(r'(\w+)\s*=\s*([^;]+);\s*(?://\s*(.*))?', content)
    evaluator = SafeExpressionEvaluator({})
    
    for name, value_expr, comment in variables_matches:
        try:
            # Atualizar evaluator com novas variáveis
            evaluator.variables = data['variables']
            
            # Avaliar expressão
            value = evaluator.evaluate(value_expr.strip())
            
            data['variables'][name] = value
            print(f"  {name} = {value} // {comment if comment else ''}")
            
        except Exception as e:
            print(f"⚠️ Erro ao processar variável {name}: {e}")

    # Atualizar evaluator final
    evaluator = SafeExpressionEvaluator(data['variables'])

    # 2. Parse Points
    print("📍 Extraindo pontos...")
    point_matches = re.findall(r'Point\((\d+)\)\s*=\s*\{([^}]+)\};', content)
    for p_id, coords_str in point_matches:
        try:
            coords = [c.strip() for c in coords_str.split(',')]
            if len(coords) >= 4:
                x = evaluator.evaluate(coords[0])
                y = evaluator.evaluate(coords[1]) 
                z = evaluator.evaluate(coords[2])
                lc = evaluator.evaluate(coords[3])
                
                point = {
                    'id': int(p_id),
                    'x': x,
                    'y': y,
                    'z': z,
                    'lc': lc
                }
                data['points'].append(point)
                print(f"  Point({p_id}): ({x:.2f}, {y:.2f}, {z:.2f}) lc={lc:.2f}")
                
        except Exception as e:
            print(f"⚠️ Erro ao processar ponto {p_id}: {e}")

    # 3. Parse Lines
    print("📏 Extraindo linhas...")
    line_matches = re.findall(r'Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\};', content)
    for l_id, p1, p2 in line_matches:
        try:
            line = {
                'id': int(l_id),
                'p1': int(p1),
                'p2': int(p2)
            }
            data['lines'].append(line)
            
        except Exception as e:
            print(f"⚠️ Erro ao processar linha {l_id}: {e}")
    
    print(f"  Extraídas {len(data['lines'])} linhas")

    # 4. Parse Curve Loops
    print("🔄 Extraindo curve loops...")
    curve_loop_matches = re.findall(r'Curve Loop\((\d+)\)\s*=\s*\{([^}]+)\};', content)
    for cl_id, lines_str in curve_loop_matches:
        try:
            lines_list = []
            for line_id_str in lines_str.split(','):
                line_id_str = line_id_str.strip()
                is_reversed = line_id_str.startswith('-')
                line_id = abs(int(line_id_str))
                lines_list.append({'id': line_id, 'reversed': is_reversed})
            
            data['curve_loops'][int(cl_id)] = lines_list
            
        except Exception as e:
            print(f"⚠️ Erro ao processar curve loop {cl_id}: {e}")
    
    print(f"  Extraídos {len(data['curve_loops'])} curve loops")

    # 5. Parse Plane Surfaces
    print("🏠 Extraindo superfícies...")
    plane_surface_matches = re.findall(r'Plane Surface\((\d+)\)\s*=\s*\{(\d+)\};', content)
    for ps_id, cl_id in plane_surface_matches:
        try:
            surface = {
                'id': int(ps_id),
                'curve_loop_id': int(cl_id)
            }
            data['plane_surfaces'].append(surface)
            
        except Exception as e:
            print(f"⚠️ Erro ao processar superfície {ps_id}: {e}")
    
    print(f"  Extraídas {len(data['plane_surfaces'])} superfícies")

    # 6. Parse Physical Groups
    print("🎯 Extraindo grupos físicos...")
    
    # Physical Surfaces
    physical_surface_matches = re.findall(r'Physical Surface\("([^"]+)",\s*(\d+)\)\s*=\s*\{([^}]+)\};', content)
    for name, p_id, surfaces_str in physical_surface_matches:
        try:
            surfaces = [int(s.strip()) for s in surfaces_str.split(',')]
            data['physical_surfaces'].append({
                'id': int(p_id),
                'name': name,
                'elements': surfaces
            })
            
        except Exception as e:
            print(f"⚠️ Erro ao processar physical surface {name}: {e}")

    # Physical Lines
    physical_line_matches = re.findall(r'Physical Line\("([^"]+)",\s*(\d+)\)\s*=\s*\{([^}]+)\};', content)
    for name, p_id, lines_str in physical_line_matches:
        try:
            lines = [int(l.strip()) for l in lines_str.split(',')]
            data['physical_lines'].append({
                'id': int(p_id),
                'name': name,
                'elements': lines
            })
            
        except Exception as e:
            print(f"⚠️ Erro ao processar physical line {name}: {e}")
    
    print(f"  Extraídos {len(data['physical_surfaces'])} physical surfaces")
    print(f"  Extraídos {len(data['physical_lines'])} physical lines")

    # 7. Parse Mesh Configurations
    print("⚙️ Extraindo configurações de malha...")
    mesh_config_matches = re.findall(r'Mesh\.(\w+)\s*=\s*([^;]+);', content)
    for param, value_str in mesh_config_matches:
        try:
            # Tentar avaliar como número
            try:
                value = evaluator.evaluate(value_str.strip())
                data['mesh_configs'][param] = value
            except:
                # Se não for número, manter como string
                data['mesh_configs'][param] = value_str.strip()
                
        except Exception as e:
            print(f"⚠️ Erro ao processar config {param}: {e}")
    
    print(f"  Extraídas {len(data['mesh_configs'])} configurações")

    return data

def create_database(db_name: str, parsed_data: Dict[str, Any]) -> None:
    """
    Cria banco de dados SQLite melhorado
    """
    print(f"🗄️ Criando banco de dados: {db_name}")
    
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Habilitar foreign keys
    cursor.execute("PRAGMA foreign_keys = ON;")

    # --- Criar Tabelas ---
    print("📋 Criando tabelas...")
    
    # Tabela de variáveis
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS variaveis (
            nome TEXT PRIMARY KEY,
            valor REAL NOT NULL,
            descricao TEXT
        );
    """)

    # Tabela de pontos
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pontos (
            id INTEGER PRIMARY KEY,
            coord_x REAL NOT NULL,
            coord_y REAL NOT NULL,
            coord_z REAL NOT NULL,
            tamanho_malha REAL
        );
    """)

    # Tabela de linhas
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS linhas (
            id INTEGER PRIMARY KEY,
            ponto_inicial_id INTEGER NOT NULL,
            ponto_final_id INTEGER NOT NULL,
            FOREIGN KEY (ponto_inicial_id) REFERENCES pontos(id),
            FOREIGN KEY (ponto_final_id) REFERENCES pontos(id)
        );
    """)

    # Tabela de superfícies
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS superficies (
            id INTEGER PRIMARY KEY,
            curve_loop_id INTEGER NOT NULL
        );
    """)

    # Tabela de curve loops
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS curve_loops (
            id INTEGER PRIMARY KEY,
            descricao TEXT
        );
    """)

    # Tabela de linhas em curve loops
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS linhas_em_curve_loops (
            curve_loop_id INTEGER NOT NULL,
            linha_id INTEGER NOT NULL,
            ordem INTEGER NOT NULL,
            invertida BOOLEAN NOT NULL DEFAULT FALSE,
            PRIMARY KEY (curve_loop_id, linha_id),
            FOREIGN KEY (curve_loop_id) REFERENCES curve_loops(id),
            FOREIGN KEY (linha_id) REFERENCES linhas(id)
        );
    """)

    # Tabela de grupos físicos
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS grupos_fisicos (
            id INTEGER PRIMARY KEY,
            nome TEXT NOT NULL UNIQUE,
            tipo TEXT NOT NULL CHECK (tipo IN ('Surface', 'Line', 'Point'))
        );
    """)

    # Tabela de elementos em grupos físicos
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS elementos_em_grupos (
            grupo_id INTEGER NOT NULL,
            elemento_id INTEGER NOT NULL,
            tipo_elemento TEXT NOT NULL,
            PRIMARY KEY (grupo_id, elemento_id),
            FOREIGN KEY (grupo_id) REFERENCES grupos_fisicos(id)
        );
    """)

    # Tabela de configurações de malha
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS configuracoes_malha (
            parametro TEXT PRIMARY KEY,
            valor TEXT NOT NULL
        );
    """)

    # --- Inserir Dados ---
    print("📊 Inserindo dados...")
    
    # Variáveis
    for nome, valor in parsed_data['variables'].items():
        cursor.execute("INSERT OR REPLACE INTO variaveis (nome, valor) VALUES (?, ?)",
                      (nome, valor))

    # Pontos
    for ponto in parsed_data['points']:
        cursor.execute("""INSERT OR REPLACE INTO pontos 
                         (id, coord_x, coord_y, coord_z, tamanho_malha) 
                         VALUES (?, ?, ?, ?, ?)""",
                      (ponto['id'], ponto['x'], ponto['y'], ponto['z'], ponto['lc']))

    # Linhas
    for linha in parsed_data['lines']:
        cursor.execute("""INSERT OR REPLACE INTO linhas 
                         (id, ponto_inicial_id, ponto_final_id) 
                         VALUES (?, ?, ?)""",
                      (linha['id'], linha['p1'], linha['p2']))

    # Curve Loops
    for cl_id, lines_info in parsed_data['curve_loops'].items():
        cursor.execute("INSERT OR REPLACE INTO curve_loops (id) VALUES (?)", (cl_id,))
        
        for i, line_data in enumerate(lines_info):
            cursor.execute("""INSERT OR REPLACE INTO linhas_em_curve_loops 
                             (curve_loop_id, linha_id, ordem, invertida) 
                             VALUES (?, ?, ?, ?)""",
                          (cl_id, line_data['id'], i + 1, line_data['reversed']))

    # Superfícies
    for superficie in parsed_data['plane_surfaces']:
        cursor.execute("INSERT OR REPLACE INTO superficies (id, curve_loop_id) VALUES (?, ?)",
                      (superficie['id'], superficie['curve_loop_id']))

    # Grupos físicos
    for grupo_type, grupos in [('Surface', parsed_data['physical_surfaces']),
                              ('Line', parsed_data['physical_lines'])]:
        for grupo in grupos:
            cursor.execute("""INSERT OR REPLACE INTO grupos_fisicos 
                             (id, nome, tipo) VALUES (?, ?, ?)""",
                          (grupo['id'], grupo['name'], grupo_type))
            
            for elemento_id in grupo['elements']:
                cursor.execute("""INSERT OR REPLACE INTO elementos_em_grupos 
                                 (grupo_id, elemento_id, tipo_elemento) 
                                 VALUES (?, ?, ?)""",
                              (grupo['id'], elemento_id, grupo_type))

    # Configurações de malha
    for param, valor in parsed_data['mesh_configs'].items():
        cursor.execute("INSERT OR REPLACE INTO configuracoes_malha (parametro, valor) VALUES (?, ?)",
                      (param, str(valor)))

    conn.commit()
    print(f"✅ Banco de dados criado com sucesso!")
    print(f"📊 Estatísticas:")
    print(f"  - {len(parsed_data['variables'])} variáveis")
    print(f"  - {len(parsed_data['points'])} pontos")  
    print(f"  - {len(parsed_data['lines'])} linhas")
    print(f"  - {len(parsed_data['curve_loops'])} curve loops")
    print(f"  - {len(parsed_data['plane_surfaces'])} superfícies")
    print(f"  - {len(parsed_data['physical_surfaces'])} physical surfaces")
    print(f"  - {len(parsed_data['physical_lines'])} physical lines")
    
    conn.close()

def query_database(db_name: str) -> None:
    """
    Executa algumas consultas de exemplo no banco
    """
    print(f"\n🔍 Executando consultas de exemplo em {db_name}...")
    
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()
    
    # Lista de pontos
    print("\n📍 Primeiros 5 pontos:")
    cursor.execute("SELECT id, coord_x, coord_y, tamanho_malha FROM pontos LIMIT 5")
    for row in cursor.fetchall():
        print(f"  Ponto {row[0]}: ({row[1]:.2f}, {row[2]:.2f}) lc={row[3]:.2f}")
    
    # Grupos físicos
    print("\n🎯 Grupos físicos:")
    cursor.execute("SELECT nome, tipo, COUNT(elemento_id) as num_elementos FROM grupos_fisicos g LEFT JOIN elementos_em_grupos e ON g.id = e.grupo_id GROUP BY g.id")
    for row in cursor.fetchall():
        print(f"  {row[0]} ({row[1]}): {row[2]} elementos")
    
    # Configurações de malha
    print("\n⚙️ Configurações de malha:")
    cursor.execute("SELECT parametro, valor FROM configuracoes_malha")
    for row in cursor.fetchall():
        print(f"  {row[0]} = {row[1]}")
    
    conn.close()

def main():
    """Função principal"""
    print("🚀 Conversor de .geo para SQLite - Versão Melhorada")
    print("=" * 50)
    
    # Usar o arquivo real da barragem
    geo_file = 'barragem1.geo'
    db_file = 'barragem_geometria.db'
    
    # Verificar se arquivo existe
    import os
    if not os.path.exists(geo_file):
        print(f"❌ Arquivo {geo_file} não encontrado!")
        return
    
    # Processar arquivo
    parsed_data = parse_geo_file(geo_file)
    
    # Criar banco de dados
    create_database(db_file, parsed_data)
    
    # Executar consultas de exemplo
    query_database(db_file)
    
    print(f"\n🎉 Processo concluído!")
    print(f"📁 Banco de dados salvo em: {db_file}")
    print(f"💡 Use um visualizador SQLite para explorar os dados")

if __name__ == "__main__":
    main() 