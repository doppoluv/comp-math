import numpy as np
import matplotlib.pyplot as plt

def df(f, x, h=1e-8):
    return (f(x + h) - f(x - h)) / (2 * h)

def separate_roots(f, a, b, n=1000):
    x = np.linspace(a, b, n + 1)
    intervals = []
    roots = []
    for i in range(n):
        if f(x[i + 1]) == 0:
            roots.append(x[i + 1]) 
        if f(x[i]) * f(x[i + 1]) < 0:
            intervals.append([x[i], x[i + 1]])
    print(intervals)
    return intervals, roots

def refine_separation(f, intervals, max_iter=1000):
    result = []
    
    for interval in intervals:
        a, b = interval
        x_test = np.linspace(a, b, 20)
        df_vals = [df(f, x) for x in x_test]
        
        if all(v > 0 for v in df_vals) or all(v < 0 for v in df_vals):
            result.append([a, b])
        else:
            mid = (a + b) / 2
            sub_intervals = []
            
            if f(a) * f(mid) < 0:
                sub_intervals.append([a, mid])
            if f(mid) * f(b) < 0:
                sub_intervals.append([mid, b])
            
            result.extend(refine_separation(f, sub_intervals, max_iter - 1))
    
    return result

def find_root_intervals(f, a=0, b=1, n=1000):
    intervals, roots = separate_roots(f, a, b, n)
    root_intervals = refine_separation(f, intervals)
    return root_intervals, roots





def bisection(f, a, b, eps=1e-8):
    iter_count = 0
    
    while (b - a) / 2 > eps:
        c = (a + b) / 2
        iter_count += 1
        
        if abs(f(c)) < eps:
            return c, iter_count
        
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    
    return (a + b) / 2, iter_count

def newton(f, a, b, eps=1e-8, max_iter=1000):
    x = (a + b) / 2
    iter_count = 0
    
    for i in range(max_iter):
        iter_count += 1
        fx = f(x)
        
        if abs(fx) < eps:
            return x, iter_count
        
        dfx = df(f, x)
        
        if abs(dfx) < 1e-12:
            return x, iter_count
        
        x_new = x - fx / dfx
        
        if abs(x_new - x) < eps:
            return x_new, iter_count
        
        x = x_new
    
    return x, iter_count

def secant(f, a, b, eps=1e-8, max_iter=1000):
    x0, x1 = a, b
    iter_count = 0
    
    for i in range(max_iter):
        iter_count += 1
        
        f0, f1 = f(x0), f(x1)
        
        if abs(f1) < eps:
            return x1, iter_count
        
        if abs(f1 - f0) < 1e-12:
            return x1, iter_count
        
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        
        if abs(x2 - x1) < eps:
            return x2, iter_count
        
        x0, x1 = x1, x2
    
    return x1, iter_count





def refine_roots(f, intervals, eps=1e-8):
    results = []
    
    for i, (a, b) in enumerate(intervals, 1):
        root_bisect, iter_bisect = bisection(f, a, b, eps)
        root_newton, iter_newton = newton(f, a, b, eps)
        root_secant, iter_secant = secant(f, a, b, eps)
        
        results.append({
            'num': i,
            'interval': (a, b),
            'bisection': {'root': root_bisect, 'iters': iter_bisect},
            'newton': {'root': root_newton, 'iters': iter_newton},
            'secant': {'root': root_secant, 'iters': iter_secant}
        })
    
    return results

def print_results(results, f):
    for res in results:
        print(f"\nКОРЕНЬ #{res['num']}:")
        print(f"Отрезок: [{res['interval'][0]}, {res['interval'][1]}]")
        print("-" * 110)
        
        print(f"{'Метод':<30} {'Корень':<30} {'f(x)':<30} {'Итерации':<10}")
        print("-" * 110)
        
        for method_name, method_key in [('Бисекция', 'bisection'), ('Ньютона', 'newton'), ('Секущих', 'secant')]:
            root = res[method_key]['root']
            iters = res[method_key]['iters']
            fx = f(root)
            print(f"{method_name:<30} {root:<30} {fx:<30.5e} {iters:<10}")

def plot_function_and_roots(f, a, b, results, title):
    x = np.linspace(a, b, 1000)
    y = np.array([f(xi) for xi in x])
    
    plt.figure(figsize=(8, 4))
    
    plt.plot(x, y, 'b-', linewidth=2, label='f(x)')
    plt.axhline(y=0, color='k', linestyle='--', linewidth=0.5, alpha=0.3)
    plt.grid(True, alpha=0.3)
    
    colors = {'bisection': 'red', 'newton': 'green', 'secant': 'orange'}
    markers = {'bisection': 'o', 'newton': 's', 'secant': '^'}
    labels = {'bisection': 'Бисекция', 'newton': 'Ньютона', 'secant': 'Секущих'}
    
    for method_key in ['bisection', 'newton', 'secant']:
        roots = [res[method_key]['root'] for res in results]
        plt.scatter(roots, [0]*len(roots), 
                    color=colors[method_key], 
                    marker=markers[method_key], 
                    s=100, 
                    label=labels[method_key],
                    zorder=5)
    
    plt.xlabel('x', fontsize=12)
    plt.ylabel('f(x)', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()






a = -7
b = 3

def f1(x):
    return 5**x - 6*x - 3

root_intervals = find_root_intervals(f1, a, b, 1000)
print(f"\nНайдено отрезков: {len(root_intervals)}")
results = refine_roots(f1, root_intervals)
print_results(results, f1)
plot_function_and_roots(f1, a, b, results, 'f(x) = 5**x - 6*x - 3')


print("\n")
print("="*110)


def f2(x):
    return x**4 - x**3 - 2*x**2 + 3*x - 3

root_intervals = find_root_intervals(f2, a, b, 1000)
print(f"\nНайдено отрезков: {len(root_intervals)}")
results = refine_roots(f2, root_intervals)
print_results(results, f2)
plot_function_and_roots(f2, a, b, results, 'f(x) = x**4 - x**3 - 2*x**2 + 3*x - 3')


print("\n")
print("="*110)


def f3(x):
    return 2*x**2 - 0.5**x - 3

root_intervals = find_root_intervals(f3, a, b, 1000)
print(f"\nНайдено отрезков: {len(root_intervals)}")
results = refine_roots(f3, root_intervals)
print_results(results, f3)
plot_function_and_roots(f3, a, b, results, 'f(x) = 2*x**2 - 0.5**x - 3')


print("\n")
print("="*110)


def f4(x):
    return x * np.log10(x + 1) - 1

root_intervals = find_root_intervals(f4, a, b, 1000)
print(f"\nНайдено отрезков: {len(root_intervals)}")
results = refine_roots(f4, root_intervals)
print_results(results, f4)
plot_function_and_roots(f4, a, b, results, 'f(x) = x * np.log10(x + 1) - 1')





l = 7.8
a = -l
b = l

def f5(x):
    return x - l*np.sin(x)

root_intervals, roots = find_root_intervals(f5, a, b, 1001)
print(f"\nНайдено отрезков: {len(root_intervals)}")
if len(roots) != 0:
    print(f"Корни на границах интервалов: {roots}")
results = refine_roots(f5, root_intervals)
print_results(results, f5)
plot_function_and_roots(f5, a, b, results, 'f(x) = x - l*np.sin(x)')