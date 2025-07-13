import sqlite3
import re

def parse_geo_file(filepath):
    """
    Parses a .geo file and extracts geometry data.
    Returns a dictionary with parsed data.
    """
    data = {
        'variables': {},
        'points': [],
        'lines': [],
        'curve_loops': {},
        'plane_surfaces': [],
        'physical_surfaces': [],
        'physical_lines': [],
        'mesh_configs': {}
    }

    with open(filepath, 'r') as f:
        content = f.read()

    # 1. Parse Variables
    # Example: lmsh = 2.5; // Tamanho característico da malha
    variables_matches = re.findall(r'(\w+)\s*=\s*([0-9.]+);\s*//\s*(.*)', content)
    for name, value, desc in variables_matches:
        try:
            data['variables'][name] = {'value': float(value), 'description': desc.strip()}
        except ValueError:
            # Handle non-numeric variables if any, though unlikely for initial block
            data['variables'][name] = {'value': value, 'description': desc.strip()}


    # 2. Parse Points
    # Example: Point(1) = {-hf, -hf, 0, lmsh2};
    point_matches = re.findall(r'Point\((\d+)\)\s*=\s*{(.*?),\s*(.*?),\s*(.*?),\s*(.*?)\};', content)
    for p_id, x, y, z, lc in point_matches:
        # Evaluate expressions for coordinates and mesh size using defined variables
        # This is a simplification; a full evaluator would be more robust.
        try:
            # Replace variable names with their values for evaluation
            # This is a basic substitution. For complex expressions, use ast.literal_eval or similar.
            x_val = eval(x.replace('hf', str(data['variables'].get('hf', {'value': 0})['value']))
                         .replace('lc1', str(data['variables'].get('lc1', {'value': 0})['value'])))
            y_val = eval(y.replace('hf', str(data['variables'].get('hf', {'value': 0})['value']))
                         .replace('hc1', str(data['variables'].get('hc1', {'value': 0})['value']))
                         .replace('hc2', str(data['variables'].get('hc2', {'value': 0})['value']))
                         .replace('hc3', str(data['variables'].get('hc3', {'value': 0})['value'])))
            z_val = eval(z) # Should be 0 in 2D
            lc_val = eval(lc.replace('lmsh', str(data['variables'].get('lmsh', {'value': 0})['value']))
                          .replace('lmsh2', str(data['variables'].get('lmsh2', {'value': 0})['value'])))

            data['points'].append({
                'id': int(p_id),
                'x': x_val,
                'y': y_val,
                'z': z_val,
                'lc': lc_val
            })
        except Exception as e:
            print(f"Warning: Could not parse point {p_id}. Error: {e}")
            print(f"X: {x}, Y: {y}, Z: {z}, LC: {lc}")

    # 3. Parse Lines
    # Example: Line(1) = {1, 2};
    line_matches = re.findall(r'Line\((\d+)\)\s*=\s*\{(\d+),\s*(\d+)\};', content)
    for l_id, p1, p2 in line_matches:
        data['lines'].append({
            'id': int(l_id),
            'p1': int(p1),
            'p2': int(p2)
        })

    # 4. Parse Curve Loops and Plane Surfaces
    # Example: Curve Loop(1) = {9, 10, 1, 11};
    # Example: Plane Surface(1) = {1};
    curve_loop_blocks = re.findall(r'Curve Loop\((\d+)\)\s*=\s*\{(.*?)\};', content)
    for cl_id, lines_str in curve_loop_blocks:
        lines_list = []
        for line_id_str in lines_str.split(','):
            line_id = int(line_id_str.strip())
            is_reversed = False
            if line_id < 0:
                is_reversed = True
                line_id = abs(line_id)
            lines_list.append({'id': line_id, 'reversed': is_reversed})
        data['curve_loops'][int(cl_id)] = lines_list

    plane_surface_matches = re.findall(r'Plane Surface\((\d+)\)\s*=\s*\{(\d+)\};', content)
    for ps_id, cl_id in plane_surface_matches:
        data['plane_surfaces'].append({
            'id': int(ps_id),
            'curve_loop_id': int(cl_id)
        })

    # 5. Parse Physical Groups
    # Example: Physical Surface("camada_material_1", 1) = {1};
    physical_surface_matches = re.findall(r'Physical Surface\("(.+?)",\s*(\d+)\)\s*=\s*\{(.*?)\};', content)
    for name, p_id, surfaces_str in physical_surface_matches:
        surfaces = [int(s.strip()) for s in surfaces_str.split(',')]
        data['physical_surfaces'].append({
            'id': int(p_id),
            'name': name,
            'elements': surfaces
        })

    # Example: Physical Line("ISOLAMENTO_PERFEITO", 11) = {10,1,2,3,4,5};
    physical_line_matches = re.findall(r'Physical Line\("(.+?)",\s*(\d+)\)\s*=\s*\{(.*?)\};', content)
    for name, p_id, lines_str in physical_line_matches:
        lines = [int(l.strip()) for l in lines_str.split(',')]
        data['physical_lines'].append({
            'id': int(p_id),
            'name': name,
            'elements': lines
        })

    # 6. Parse Mesh Configurations
    # Example: Mesh.MeshSizeMin = lmsh;
    mesh_config_matches = re.findall(r'Mesh\.(\w+)\s*=\s*(.*?);', content)
    for param, value_str in mesh_config_matches:
        # Try to resolve variable values for mesh configs
        resolved_value = value_str.strip()
        if resolved_value in data['variables']:
            resolved_value = str(data['variables'][resolved_value]['value'])
        data['mesh_configs'][param] = resolved_value


    return data

def create_and_populate_db(db_name, parsed_data):
    """
    Creates an SQLite database and populates it with parsed .geo data.
    """
    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    # Disable foreign key checks for easier bulk insertion
    cursor.execute("PRAGMA foreign_keys = OFF;")

    # --- Create Tables ---
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Variaveis (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nome TEXT NOT NULL UNIQUE,
            valor REAL NOT NULL,
            descricao TEXT
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Pontos (
            id INTEGER PRIMARY KEY,
            coord_x REAL NOT NULL,
            coord_y REAL NOT NULL,
            coord_z REAL NOT NULL,
            tamanho_malha REAL
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Linhas (
            id INTEGER PRIMARY KEY,
            ponto_inicial_id INTEGER NOT NULL,
            ponto_final_id INTEGER NOT NULL,
            FOREIGN KEY (ponto_inicial_id) REFERENCES Pontos(id),
            FOREIGN KEY (ponto_final_id) REFERENCES Pontos(id)
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Superficies (
            id INTEGER PRIMARY KEY,
            id_curve_loop INTEGER NOT NULL UNIQUE
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS GruposFisicos (
            id INTEGER PRIMARY KEY,
            nome TEXT NOT NULL UNIQUE,
            tipo TEXT NOT NULL
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS ConfiguracoesMalha (
            parametro TEXT PRIMARY KEY,
            valor TEXT NOT NULL
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Linhas_em_CurveLoops (
            id_curve_loop INTEGER NOT NULL,
            id_linha INTEGER NOT NULL,
            ordem INTEGER NOT NULL,
            sentido_inverso BOOLEAN NOT NULL DEFAULT FALSE,
            PRIMARY KEY (id_curve_loop, id_linha),
            FOREIGN KEY (id_curve_loop) REFERENCES Superficies(id),
            FOREIGN KEY (id_linha) REFERENCES Linhas(id)
        );
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Elementos_em_GruposFisicos (
            id_grupo_fisico INTEGER NOT NULL,
            id_elemento INTEGER NOT NULL,
            tipo_elemento TEXT NOT NULL,
            PRIMARY KEY (id_grupo_fisico, id_elemento, tipo_elemento),
            FOREIGN KEY (id_grupo_fisico) REFERENCES GruposFisicos(id)
        );
    """)

    # --- Populate Tables ---

    # Variables
    for name, info in parsed_data['variables'].items():
        cursor.execute("INSERT INTO Variaveis (nome, valor, descricao) VALUES (?, ?, ?)",
                       (name, info['value'], info['description']))

    # Points
    for p in parsed_data['points']:
        cursor.execute("INSERT INTO Pontos (id, coord_x, coord_y, coord_z, tamanho_malha) VALUES (?, ?, ?, ?, ?)",
                       (p['id'], p['x'], p['y'], p['z'], p['lc']))

    # Lines
    for l in parsed_data['lines']:
        cursor.execute("INSERT INTO Linhas (id, ponto_inicial_id, ponto_final_id) VALUES (?, ?, ?)",
                       (l['id'], l['p1'], l['p2']))

    # Surfaces
    for s in parsed_data['plane_surfaces']:
        cursor.execute("INSERT INTO Superficies (id, id_curve_loop) VALUES (?, ?)",
                       (s['id'], s['curve_loop_id']))

    # Linhas_em_CurveLoops
    for curve_loop_id, lines_info in parsed_data['curve_loops'].items():
        for i, line_data in enumerate(lines_info):
            cursor.execute("INSERT INTO Linhas_em_CurveLoops (id_curve_loop, id_linha, ordem, sentido_inverso) VALUES (?, ?, ?, ?)",
                           (curve_loop_id, line_data['id'], i + 1, line_data['reversed']))

    # Physical Groups
    for pg_type, pg_list in [('Surface', parsed_data['physical_surfaces']),
                             ('Line', parsed_data['physical_lines'])]:
        for pg in pg_list:
            cursor.execute("INSERT INTO GruposFisicos (id, nome, tipo) VALUES (?, ?, ?)",
                           (pg['id'], pg['name'], pg_type))
            for element_id in pg['elements']:
                cursor.execute("INSERT INTO Elementos_em_GruposFisicos (id_grupo_fisico, id_elemento, tipo_elemento) VALUES (?, ?, ?)",
                               (pg['id'], element_id, pg_type))

    # Mesh Configurations
    for param, value in parsed_data['mesh_configs'].items():
        cursor.execute("INSERT INTO ConfiguracoesMalha (parametro, valor) VALUES (?, ?)",
                       (param, value))

    # Enable foreign key checks again
    cursor.execute("PRAGMA foreign_keys = ON;")

    conn.commit()
    conn.close()
    print(f"Database '{db_name}' created and populated successfully.")

if __name__ == "__main__":
    # Example usage:
    geo_file = 'meu_arquivo.geo' # <-- Mude para o nome do seu arquivo .geo
    db_file = 'geometria.db'     # <-- Nome do arquivo de banco de dados SQLite a ser criado

    # **PASTE YOUR .geo CONTENT HERE TO CREATE 'meu_arquivo.geo'**
    # For testing, you can paste the .geo content you provided earlier into a file named 'meu_arquivo.geo'
    # Or, replace this block with direct file operations.
    geo_content_example = """
// Definição de pontos
lmsh = 2.5; // Tamanho característico da malha
lmsh2 = 3.0; // Tamanho característico da malha
hf = 3.0; // profundidade da fundação
hc1 = 2.5; // altura da camada 1
hc2 = 5.0; // altura da camada 2
hc3 = 7.5; // altura da camada 3
lc1 = 5.0; // largura da camada 1
lc2 = 4.0; // largura da camada 2
lc3 = 3.0; // largura da camada 3
lcf = 1.0; // largura do concreto de face
lf = lc1 + hf; // largura da fundação

// Coordenadas dos pontos
Point(1) = {-hf, -hf, 0, lmsh2};
Point(2) = {0, -hf, 0, lmsh2};
Point(3) = {lcf, -hf, 0, lmsh2};
Point(4) = {lc1, -hf, 0, lmsh2};
Point(5) = {lf, -hf, 0, lmsh2};
Point(6) = {lf, 0, 0, lmsh2};
Point(7) = {lc1, 0, 0, lmsh};
Point(8) = {lcf, 0, 0, lmsh};
Point(9) = {0, 0, 0, lmsh};
Point(10) = {-hf, 0, 0, lmsh2};
Point(11) = {0, hc1, 0, lmsh};
Point(12) = {lcf, hc1, 0, lmsh};
Point(13) = {lc2, hc1, 0, lmsh};
Point(14) = {lc3, hc2, 0, lmsh};
Point(15) = {lcf, hc2, 0, lmsh};
Point(16) = {0, hc2, 0, lmsh};
Point(17) = {0, hc3, 0, lmsh};
Point(18) = {lcf, hc3, 0, lmsh};
Point(19) = {lc3, hc3, 0, lmsh};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};
Line(11) = {2, 9};
Line(12) = {8, 3};
Line(13) = {4, 7};
Line(14) = {7, 13};
Line(15) = {13, 12};
Line(16) = {12, 11};
Line(17) = {11, 9};
Line(18) = {8, 12};
Line(19) = {11, 16};
Line(20) = {16, 15};
Line(21) = {15, 14};
Line(22) = {14, 13};
Line(23) = {12, 15};
Line(24) = {14, 19};
Line(25) = {19, 18};
Line(26) = {18, 17};
Line(27) = {17, 16};
Line(28) = {15, 18};

Curve Loop(1) = {9, 10, 1, 11};
Plane Surface(1) = {1};
Curve Loop(2) = {8, -11, 2, -12};
Plane Surface(2) = {2};
Curve Loop(3) = {7, 12, 3, 13};
Plane Surface(3) = {3};
Curve Loop(4) = {6, -13, 4, 5};
Plane Surface(4) = {4};
Curve Loop(5) = {7, 18, -15, -14};
Plane Surface(5) = {5};
Curve Loop(6) = {16, 17, -8, 18};
Plane Surface(6) = {6};
Curve Loop(7) = {23, -20, -19, -16};
Plane Surface(7) = {7};
Curve Loop(8) = {21, 22, 15, 23};
Plane Surface(8) = {8};
Curve Loop(9) = {24, 25, -28, 21};
Plane Surface(9) = {9};
Curve Loop(10) = {26, 27, 20, 28};
Plane Surface(10) = {10};

// Physical Groups para as camadas de material
Physical Surface("camada_material_1", 1) = {1};
Physical Surface("camada_material_2", 2) = {2};
Physical Surface("camada_material_3", 3) = {3};
Physical Surface("camada_material_4", 4) = {4};
Physical Surface("camada_material_5", 5) = {6};
Physical Surface("camada_material_6", 6) = {5};
Physical Surface("camada_material_7", 7) = {7};
Physical Surface("camada_material_8", 8) = {8};
Physical Surface("camada_material_9", 9) = {10};
Physical Surface("camada_material_10", 10) = {9};

Physical Line("ISOLAMENTO_PERFEITO", 11) = {10,1,2,3,4,5};
Physical Line("FUNDACAO_TOPO", 12) = {9,6};
Physical Line("FACE_MONTANTE_1", 13) = {17};
Physical Line("FACE_JUSANTE_1", 14) = {14};
Physical Line("FACE_MONTANTE_2", 15) = {19};
Physical Line("FACE_JUSANTE_2", 16) = {22};
Physical Line("FACE_MONTANTE_3", 17) = {27};
Physical Line("FACE_JUSANTE_3", 18) = {24};
Physical Line("FACE_TOPO", 19) = {26,25};
Physical Line("interface_1_2", 20) = {16,15};
Physical Line("interface_2_3", 21) = {20,21};

// Configurações da malha
// Definindo Transfinite nas linhas (mínimo possível)
Transfinite Line {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28} = 2;

// Aplicando Transfinite nas superfícies
Transfinite Surface {1,2,3,4,5,6,7,8,9,10};

// Recombinando as superfícies para quadriláteros perfeitos
Recombine Surface {1,2,3,4,5,6,7,8,9,10};

// Mudando o algoritmo para estruturado
Mesh.MeshSizeMin = lmsh;
Mesh.MeshSizeMax = lmsh2;
Mesh.RecombineAll = 1;
Mesh.Algorithm = 1;

Mesh.MshFileVersion = 2.2;
    """
    with open(geo_file, 'w') as f:
        f.write(geo_content_example)
    # --- END PASTE BLOCK ---

    print(f"Parsing '{geo_file}'...")
    parsed_data = parse_geo_file(geo_file)

    print(f"Creating and populating database '{db_file}'...")
    create_and_populate_db(db_file, parsed_data)

    print("\nProcess finished.")
    print(f"You can now explore '{db_file}' using a SQLite browser or Python.")