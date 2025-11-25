import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def f(x):
    return np.sqrt(x**2)

def generate_data(a, b, N):
    # # Равноммерные узлы
    # h = (b - a) / N
    # xi = np.array([a + i * h for i in range(N + 1)])
    # fi = f(xi)

    # Узлы Чебышёва (для Лагранжа на отрезке [-1;1])
    i = np.arange(N + 1)
    xi = np.cos((2 * i + 1) * np.pi / (2 * (N + 1)))
    fi = f(xi)

    return xi, fi

def lagrange_interpolation(x, xi, fi):
    n = len(xi)
    result = 0.0
    for i in range(n):
        term = fi[i]
        for j in range(n):
            if j != i:
                term *= (x - xi[j]) / (xi[i] - xi[j])
        result += term
    return result

def cubic_spline(xi, fi):
    n = len(xi) - 1
    h = np.diff(xi)
    alpha = np.zeros(n)
    c = np.zeros(n + 1)
    d = np.zeros(n + 1)
    b = np.zeros(n)

    for i in range(1, n):
        alpha[i] = (6 / h[i] * (fi[i + 1] - fi[i]) - 6 / h[i - 1] * (fi[i] - fi[i - 1]))

    l = np.zeros(n + 1)
    mu = np.zeros(n + 1)
    z = np.zeros(n + 1)
    l[0] = 1
    mu[0] = 0
    z[0] = 0

    for i in range(1, n):
        l[i] = 2 * (xi[i + 1] - xi[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[n] = 1
    z[n] = 0
    c[n] = 0

    for i in range(n - 1, -1, -1):
        c[i] = z[i] - mu[i] * c[i + 1]
        b[i] = (fi[i + 1] - fi[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 6
        d[i] = (c[i + 1] - c[i]) / (h[i])

    a = fi[:-1]
    return a, b, c[:-1], d

def eval_spline(x, xi, a, b, c, d):
    n = len(xi) - 1
    for i in range(n):
        if xi[i] <= x <= xi[i + 1]:
            dx = x - xi[i]
            return a[i] + b[i] * dx + c[i] * dx**2 / 2 + d[i] * dx**3 / 6
    return None

def plot_interpolation(a, b, N, num_points=1000):
    xi, fi = generate_data(a, b, N)

    print(f"Таблица для N={N}:")
    for i in range(len(xi)):
        print(f"xi[{i}] = {xi[i]:.4f}, fi[{i}] = {fi[i]:.4f}")
    print("\n")

    x_eval = np.linspace(a, b, num_points)

    y_original = f(x_eval)

    y_lagrange = np.array([lagrange_interpolation(x, xi, fi) for x in x_eval])

    # a, b, c, d = cubic_spline(xi, fi)
    # y_spline_manual = np.array([eval_spline(x, xi, a, b, c, d) for x in x_eval])

    # cs = CubicSpline(xi, fi, bc_type='natural')
    # y_spline_lib = cs(x_eval)



    err_l = np.abs(y_lagrange - y_original)
    # err_spl = np.abs(y_spline_manual - y_original)
    print(f"Ошибка Ланганж: {np.max(err_l)}")
    # print(f"Ошибка сплайна: {np.max(err_spl)}")
    print("\n")



    plt.figure(figsize=(10, 6))
    plt.plot(x_eval, y_original, 'k-', label='f(x)')
    plt.plot(x_eval, y_lagrange, 'b--', label='Интерполяционный полином')
    # plt.plot(x_eval, y_spline_manual, 'r-.', label='Кубический сплайн')
    # plt.plot(x_eval, y_spline_lib, 'g--', label='Библиотечный сплайн')
    plt.plot(xi, fi, 'yo', label='Узлы')
    plt.title(f'N={N}')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.show()

a = -1.0
b = 1.0
Ns = [2, 5, 10, 20]
for N in Ns:
    plot_interpolation(a, b, N)

