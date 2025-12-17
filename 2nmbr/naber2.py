import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Исходные данные
# ===============================
f = 6.5*10**9                 # частота, Гц
c = 3*10**8                   # скорость света, м/с
lambd = c / f               # длина волны, м
k = 2 * np.pi / lambd       # волновое число

# По условию: 2l / lambda = 0.01
l = 0.005 * lambd           # длина одного плеча вибратора

# ===============================
# Угловая сетка
# ===============================
theta = np.linspace(1e-6, np.pi - 1e-6, 2000)

# ===============================
# Напряженность поля (формула 1)
# ===============================
E_theta = (np.cos(k * l * np.cos(theta)) - np.cos(k * l)) / np.sin(theta)

# ===============================
# Нормированная ДН по полю (формула 3)
# ===============================
F = np.abs(E_theta)
F = F / np.max(F)

# ===============================
# Расчет максимального КНД (формула 2)
# ===============================
integral = np.trapz(F**2 * np.sin(theta), theta)
D_max = 4 * np.pi / (2 * np.pi * integral)
D_max_db = 10 * np.log10(D_max)

print("Максимальный КНД:")
print(f"D_max = {D_max:.4f}")
print(f"D_max = {D_max_db:.2f} дБ")

# ===============================
# КНД как функция угла (формула 4)
# ===============================
D_theta = D_max * F**2
D_theta_db = 10 * np.log10(D_theta)

# ===============================
# Графики в декартовой системе
# ===============================
plt.figure()
plt.plot(theta * 180 / np.pi, D_theta)
plt.xlabel("θ, град")
plt.ylabel("D(θ)")
plt.title("КНД в разах")
plt.grid()
plt.show()

plt.figure()
plt.plot(theta * 180 / np.pi, D_theta_db)
plt.xlabel("θ, град")
plt.ylabel("D(θ), дБ")
plt.title("КНД в дБ")
plt.grid()
plt.show()

# ===============================
# Полярные диаграммы
# ===============================
plt.figure()
ax = plt.subplot(111, polar=True)
ax.plot(theta, D_theta)
ax.set_title("Полярная диаграмма КНД (разы)")
plt.show()

plt.figure()
ax = plt.subplot(111, polar=True)
ax.plot(theta, D_theta_db)
ax.set_title("Полярная диаграмма КНД (дБ)")
plt.show()