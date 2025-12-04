import numpy as np
import matplotlib.pyplot as plt

def exact_solution(x, y0):
    return y0 * np.exp(-x)


def euler(f, y0, x_grid):
    y = np.zeros_like(x_grid)
    y[0] = y0
    tau = x_grid[1] - x_grid[0]
    for n in range(len(x_grid)-1):
        y[n+1] = y[n] + tau * f(y[n])
    return y

def rk2(f, y0, x_grid):
    y = np.zeros_like(x_grid)
    y[0] = y0
    tau = x_grid[1] - x_grid[0]
    for n in range(len(x_grid)-1):
        k1 = f(y[n])
        k2 = f(y[n] + tau/2 * k1)
        y[n+1] = y[n] + tau * k2
    return y

def rk4(f, y0, x_grid):
    y = np.zeros_like(x_grid)
    y[0] = y0
    tau = x_grid[1] - x_grid[0]
    for n in range(len(x_grid)-1):
        yn = y[n]
        k1 = f(yn)
        k2 = f(yn + tau/2*k1)
        k3 = f(yn + tau/2*k2)
        k4 = f(yn + tau*k3)
        y[n+1] = yn + tau/6*(k1 + 2*k2 + 2*k3 + k4)
    return y



def runge_error(y_coarse, y_fine, p):
    return abs(y_coarse - y_fine) / (2**p - 1)

def empirical_order(y_h, y_h2, y_h4):
    num = abs(y_h - y_h2)
    den = abs(y_h2 - y_h4)
    return np.log2(num/den) if den > 1e-16 else float('nan')



def check_stability_by_lax(y0=-1.0, T=5.0, N=120):
    x = np.linspace(0, T, N+1)
    tau = x[1] - x[0]

    lambda_euler = 1 + (-tau)
    abs_euler = abs(lambda_euler)

    lambda_rk2 = 1 + (-tau) + (-tau)**2/2
    abs_rk2 = abs(lambda_rk2)

    lambda_rk4 = 1 + (-tau) + (-tau)**2/2 + (-tau)**3/6 + (-tau)**4/24
    abs_rk4 = abs(lambda_rk4)

    print(f"{'Метод':<22} {'|lambda|':>12} {'Устойчив'}")
    print("-"*85)
    print(f"{'Эйлер':<22} {abs_euler:12.4f} {'Да' if abs_euler <= 1 else 'Нет'}")
    print(f"{'Рунге-Кутта 2':<22} {abs_rk2:12.4f} {'Да' if abs_rk2 <= 1 else 'Нет'}")
    print(f"{'Рунге-Кутта 4':<22} {abs_rk4:12.4f} {'Да' if abs_rk4 <= 1 else 'Нет'}")



def numeric_order(y0=-1.0, T=5.0, N0=50):
    print("="*80)
    f = lambda y: -y
    Ns = [N0, 2*N0, 4*N0]
    y_exact_T = exact_solution(T, y0)

    methods = [
        ("Эйлер", euler, 1),
        ("Рунге-Кутта 2", rk2, 2),
        ("Рунге-Кутта 4", rk4, 4)
    ]

    print(f"{'Метод':<22} {'N':>6} {'tao':>10} {'y(T)':>16} {'err_exact':>12} {'err_runge':>12} {'p_факт':>8}")
    print("-"*80)

    for name, func, p_th in methods:
        y_vals = []
        for N in Ns:
            x = np.linspace(0, T, N+1)
            y = func(f, y0, x)
            y_vals.append(y[-1])


        y_h, y_h2, y_h4 = y_vals
        err_exact = abs(y_h2 - y_exact_T)
        err_runge = runge_error(y_h, y_h2, p_th)
        p_emp = empirical_order(y_h, y_h2, y_h4)

        print(f"{name:<22} {Ns[1]:6} {T/Ns[1]:10.6f} {y_h2:16.10f}"
              f"{err_exact:12.2e} {err_runge:12.2e} {p_emp:8.3f}")



def plot_schemes(y0=-1.0, T=6.0, N=120):
    x = np.linspace(0, T, N+1)
    h = x[1]-x[0]
    f = lambda y: -y

    y_ex = exact_solution(x, y0)
    y_eu = euler(f, y0, x)
    y_tr = rk2(f, y0, x)
    y_r4 = rk4(f, y0, x)

    plt.figure(figsize=(11,6))
    plt.plot(x, y_ex, 'k-',  lw=2.5, label='Точное решение y(x)=y_0 * e^(-x)')
    plt.plot(x, y_eu, 'r--', lw=2,   label='Эйлер')
    plt.plot(x, y_tr, 'b-.', lw=2,   label='Рунге-Кутта 2')
    plt.plot(x, y_r4, 'g:',  lw=3,   label='Рунге-Кутта 4')
    plt.title(f'y` = -y,  y(0)={y0},  N={N},  tao={h:.5f}')
    plt.xlabel('x'); plt.ylabel('y(x)')
    plt.legend(); plt.grid(alpha=0.4)
    plt.tight_layout(); plt.show()

    plt.figure(figsize=(11,6))
    plt.semilogy(x, abs(y_eu-y_ex), 'r--', lw=2, label='Эйлер')
    plt.semilogy(x, abs(y_tr-y_ex), 'b-.', lw=2, label='RK2')
    plt.semilogy(x, abs(y_r4-y_ex), 'g:',  lw=3, label='RK4')
    plt.title('Погрешности решений')
    plt.xlabel('x'); plt.ylabel('|y_n - y(x_n)|')
    plt.legend(); plt.grid(True,which='both',alpha=0.4)
    plt.tight_layout(); plt.show()

def runge_analysis_along_grid(y0=-1.0, T=6.0, N_coarse=60):
    f = lambda y: -y
    Ns = [N_coarse, 2*N_coarse, 4*N_coarse]
    x_fine = np.linspace(0, T, Ns[2]+1)

    methods = [
        ("Эйлера", euler, 1, 'r--'),
        ("RK2", rk2, 2, 'b-.'),
        ("RK4", rk4, 4, 'g-')
    ]

    plt.figure(figsize=(13, 9))

    for name, method, p_th, style in methods:
        x_h  = np.linspace(0, T, Ns[0]+1)
        x_h2 = np.linspace(0, T, Ns[1]+1)
        x_h4 = x_fine

        y_h  = method(f, y0, x_h)
        y_h2 = method(f, y0, x_h2)
        y_h4 = method(f, y0, x_h4)

        y_coarse = y_h
        y_mid    = y_h2[::2]
        y_fine   = y_h4[::4]

        err_runge = np.abs(y_coarse - y_mid) / (2**p_th - 1)

        num = np.abs(y_coarse - y_mid)
        den = np.abs(y_mid - y_fine)
        p_observed = np.full_like(err_runge, np.nan)
        mask = den > 1e-15
        p_observed[mask] = np.log2(num[mask] / den[mask])

        plt.subplot(2, 1, 1)
        plt.semilogy(x_h, err_runge, style, linewidth=2.2, label=f'{name}')

        plt.subplot(2, 1, 2)
        plt.plot(x_h, p_observed, style, linewidth=2.2, label=f'{name}')

    plt.subplot(2, 1, 1)
    plt.title('Погрешности в каждом узле')
    plt.ylabel('Оценка')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=11)

    plt.subplot(2, 1, 2)
    plt.title('Порядок в каждом узле')
    plt.xlabel('$x$')
    plt.ylabel('$p(x_n)$')
    plt.ylim(0, 5)
    plt.axhline(1, color='gray', linestyle=':', alpha=0.7)
    plt.axhline(2, color='gray', linestyle=':', alpha=0.7)
    plt.axhline(4, color='gray', linestyle=':', alpha=0.7)
    plt.grid(True, ls="--", alpha=0.5)
    plt.legend(fontsize=11)

    plt.tight_layout()
    plt.show()




check_stability_by_lax(y0=-1.0, T=5.0, N=120)
numeric_order(y0=-1.0, T=5.0, N0=60)
plot_schemes(y0=-1.0, T=5.0, N=120)

runge_analysis_along_grid(y0=-1.0, T=5.0, N_coarse=60)