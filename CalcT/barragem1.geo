// Definição de pontos
lmsh = 2.0; // Tamanho característico da malha
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
