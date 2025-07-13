#!/usr/bin/env python3
"""
Problema Exotérmico 2D - Bloco 1m x 1m (Isolamento Perfeito)
Equações acopladas: temperatura (Tp) e geração de calor (teq)
Condições de contorno: isolamento perfeito (Neumann homogênea)
Tempo de simulação: 10 dias
"""

from fenics import *
import numpy as np
import matplotlib.pyplot as plt
import os

# Parâmetros do problema (do arquivo .pde)
Tad_inf1 = 16.40021368590
a1 = 0.42998953984 * 24 * 3600  # Convertendo para segundos
c1 = 1.93329591299
Tad_inf2 = 13.59978631410
a2 = 1.34507960663 * 24 * 3600  # Convertendo para segundos
c2 = 1.57098129317

# Propriedades materiais (valores típicos)
k = 2.6  # Condutividade térmica (W/m·K)
rho = 2400  # Densidade (kg/m³)
ce = 900  # Calor específico (J/kg·K)
EaR = 4000  # Energia de ativação/R (K)

# Parâmetros de tempo
T_total = 10 * 24  # 10 dias em horas
T_split = 48       # 2 dias em horas

dt1 = 2.4          # dt pequeno (h) até 48h
n1 = int(T_split / dt1)
dt2 = 12.0          # dt maior (h) após 48h
n2 = int((T_total - T_split) / dt2)

print(f"Simulação: {T_total} horas (10 dias)")
print(f"Fase 1: dt={dt1}h, {n1} passos até {T_split}h")
print(f"Fase 2: dt={dt2}h, {n2} passos até {T_total}h")

# Criar malha
def mesh_and_spaces():
    mesh = RectangleMesh(Point(0, 0), Point(1, 1), 20, 20)
    V = FunctionSpace(mesh, 'P', 1)
    elemento_misto = MixedElement([V.ufl_element(), V.ufl_element()])
    W = FunctionSpace(mesh, elemento_misto)
    return mesh, V, W

mesh, V, W = mesh_and_spaces()
print(f"Malha criada: {mesh.num_vertices()} vértices, {mesh.num_cells()} células")

# Condições de contorno: isolamento perfeito (sem condições de contorno)
# Para isolamento perfeito, não aplicamos condições de contorno Dirichlet

# Funções para solução atual e anterior
u = Function(W)  # (Tp, teq)
u_n = Function(W)  # (Tp_n, teq_n)

# Extrair componentes
Tp, teq = split(u)
Tp_n, teq_n = split(u_n)

# Funções de teste
v_Tp, v_teq = TestFunctions(W)

# Condições iniciais
Tp_0 = Constant(25.0)  # Temperatura inicial: 25°C
teq_0 = Constant(0.0)  # Tempo equivalente inicial: 0
assigner = FunctionAssigner(W, [V, V])
Tp_init = interpolate(Tp_0, V)
teq_init = interpolate(teq_0, V)
assigner.assign(u_n, [Tp_init, teq_init])

# Definição de Q
def Q_function(teq_val):
    """Calcula Q baseado no valor de teq"""
    return rho * ce * (Tad_inf1 * teq_val**c1 / (a1**c1 + teq_val**c1) + 
                      Tad_inf2 * teq_val**c2 / (a2**c2 + teq_val**c2))

# Derivada de Q em relação a teq
def dQ_dteq(teq_val):
    """Derivada de Q em relação a teq"""
    term1 = Tad_inf1 * c1 * teq_val**(c1-1) * a1**c1 / (a1**c1 + teq_val**c1)**2
    term2 = Tad_inf2 * c2 * teq_val**(c2-1) * a2**c2 / (a2**c2 + teq_val**c2)**2
    return rho * ce * (term1 + term2)

# Fator de Arrhenius
def arrhenius_factor(Tp_val):
    """Fator de Arrhenius"""
    return exp(EaR * (1/298.15 - 1/(Tp_val + 273.15)))

# Resultados
os.makedirs('resultados_isolamento', exist_ok=True)
times = []
Tp_center = []
teq_center = []
Q_values = []
center_point = Point(0.5, 0.5)
time = 0.0
step = 0
print("\nIniciando simulação...")

# Fase 1: dt pequeno
for n in range(n1):
    dt = dt1
    F_Tp = (rho * ce * (Tp - Tp_n) / (dt*3600) * v_Tp * dx + 
            k * dot(grad(Tp), grad(v_Tp)) * dx - 
            dQ_dteq(teq_n) * (teq - teq_n) / (dt*3600) * v_Tp * dx)
    F_teq = ((teq - teq_n) / (dt*3600) * v_teq * dx - 
             arrhenius_factor(Tp_n) * v_teq * dx)
    F = F_Tp + F_teq
    solve(F == 0, u)
    u_n.assign(u)
    time += dt
    step += 1
    if step % 1 == 0:
        Tp_val = u.sub(0)(center_point)
        teq_val = u.sub(1)(center_point)
        Q_val = Q_function(teq_val)
        times.append(time)
        Tp_center.append(Tp_val)
        teq_center.append(teq_val)
        Q_values.append(Q_val)
        print(f"Tempo: {time:.1f} h, Tp: {Tp_val:.2f}°C, teq: {teq_val/3600:.2f}h, Q: {Q_val:.2e} W/m³")
        Tp_file = File(f"resultados_isolamento/Tp_{step:03d}.pvd")
        teq_file = File(f"resultados_isolamento/teq_{step:03d}.pvd")
        Tp_file << u.sub(0)
        teq_file << u.sub(1)

# Fase 2: dt maior
for n in range(n2):
    dt = dt2
    F_Tp = (rho * ce * (Tp - Tp_n) / (dt*3600) * v_Tp * dx + 
            k * dot(grad(Tp), grad(v_Tp)) * dx - 
            dQ_dteq(teq_n) * (teq - teq_n) / (dt*3600) * v_Tp * dx)
    F_teq = ((teq - teq_n) / (dt*3600) * v_teq * dx - 
             arrhenius_factor(Tp_n) * v_teq * dx)
    F = F_Tp + F_teq
    solve(F == 0, u)
    u_n.assign(u)
    time += dt
    step += 1
    if step % 1 == 0:
        Tp_val = u.sub(0)(center_point)
        teq_val = u.sub(1)(center_point)
        Q_val = Q_function(teq_val)
        times.append(time)
        Tp_center.append(Tp_val)
        teq_center.append(teq_val)
        Q_values.append(Q_val)
        print(f"Tempo: {time:.1f} h, Tp: {Tp_val:.2f}°C, teq: {teq_val/3600:.2f}h, Q: {Q_val:.2e} W/m³")
        Tp_file = File(f"resultados_isolamento/Tp_{step:03d}.pvd")
        teq_file = File(f"resultados_isolamento/teq_{step:03d}.pvd")
        Tp_file << u.sub(0)
        teq_file << u.sub(1)

print("\nSimulação concluída!")

# Plotar resultados
plt.figure(figsize=(15, 5))

# Temperatura
plt.subplot(1, 3, 1)
plt.plot(times, Tp_center, 'b-', linewidth=2)
plt.xlabel('Tempo (horas)')
plt.ylabel('Temperatura (°C)')
plt.title('Temperatura no Centro')
plt.grid(True)

# Tempo equivalente
plt.subplot(1, 3, 2)
plt.plot(times, np.array(teq_center)/3600, 'r-', linewidth=2)
plt.xlabel('Tempo (horas)')
plt.ylabel('Tempo Equivalente (horas)')
plt.title('Tempo Equivalente no Centro')
plt.grid(True)

# Geração de calor
plt.subplot(1, 3, 3)
plt.plot(times, Q_values, 'g-', linewidth=2)
plt.xlabel('Tempo (horas)')
plt.ylabel('Q (W/m³)')
plt.title('Geração de Calor no Centro')
plt.grid(True)

plt.tight_layout()
plt.savefig('resultados_isolamento/resultados_exotermico_isolamento.png', dpi=300, bbox_inches='tight')
plt.show()

# Salvar dados numéricos
np.savetxt('resultados_isolamento/dados_centro.txt', 
           np.column_stack([times, Tp_center, teq_center, Q_values]),
           header='Tempo(h) Temperatura(C) TempoEquivalente(s) Q(W/m3)',
           fmt='%.6e')

print("\n✅ Simulação exotérmica (Isolamento Perfeito) concluída com sucesso!")
print("📊 Gráficos salvos em 'resultados_isolamento/resultados_exotermico_isolamento.png'")
print("📁 Arquivos de saída salvos na pasta 'resultados_isolamento/'") 