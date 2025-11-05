import numpy as np

def cubic_function(x, a, b, c, d):
    return a * x ** 3 + b * x ** 2 + c * x + d

def find_root_bounds(a, b, c, d):
    A = max(abs(b), abs(c), abs(d))
    R = 1 + A / abs(a) if a != 0 else 0

    B = max(abs(a), abs(b), abs(c))
    r = 1 / (1 + B / abs(d)) if d != 0 else 0

    return r, R

def separate_roots(f, a, b, n=10000):
    x = np.linspace(a, b, n + 1)
    intervals = []

    for i in range(n):
        if f(x[i]) * f(x[i + 1]) <= 0:
            intervals.append([x[i], x[i + 1]])

    return intervals

def bisection_method(f, a, b, eps=1e-8):
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



def analyze_cubic_equation(a, b, c, d):
    print(f"\n{a}*x**3 + {b}*x**2 + {c}*x + {d} = 0")

    r, R = find_root_bounds(a, b, c, d)
    print(f"\n1. Аналитическое отделение корней:")
    print(f"   Все корни находятся в кольце с радиусами: r = {r:.6f}, R = {R:.6f}")
    print(f"   Интервалы поиска: [-{R:.6f}, -{r:.6f}], [{r:.6f}, {R:.6f}]")


    def f(x):
        return cubic_function(x, a, b, c, d)

    search_interval1 = [-R, -r]
    intervals1 = separate_roots(f, search_interval1[0], search_interval1[1])
    search_interval2 = [r, R]
    intervals2 = separate_roots(f, search_interval2[0], search_interval2[1])
    intervals = intervals1 + intervals2

    print(f"\n2. Отделение корней методом бисекции:")
    print(f"   Найдено интервалов с корнями: {len(intervals)}")

    for i, interval in enumerate(intervals, 1):
        print(f"   Интервал {i}: [{interval[0]:.6f}, {interval[1]:.6f}]")

    print(f"\n3. Уточнение корней методом бисекции:")
    print("-" * 70)
    print(f"{'№':<3} {'Корень':<25} {'f(корень)':<25} {'Итерации':<10}")
    print("-" * 70)

    roots = []
    for i, interval in enumerate(intervals, 1):
        root, iterations = bisection_method(f, interval[0], interval[1])
        roots.append(root)
        print(f"{i:<3} {root:<25.8f} {f(root):<25.2e} {iterations:<10}")




a3, b3, c3, d3 = 2, -3, -1, 1
analyze_cubic_equation(a3, b3, c3, d3)