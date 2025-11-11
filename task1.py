import numpy as np

def cubic_function(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

def find_root_bounds(a, b, c, d):
    if a == 0:
        raise ValueError("a не должен быть равен нулю")

    A = max(abs(b), abs(c), abs(d))
    R = 1 + A / abs(a)

    B = max(abs(a), abs(b), abs(c))
    r = 1 / (1 + B / abs(d)) if d != 0 else 0

    return r, R


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
    return intervals, roots


def refine_separation(f, intervals, max_iter=10000):
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

            if max_iter > 0:
                result.extend(refine_separation(f, sub_intervals, max_iter - 1))
            else:
                result.append([a, b])

    return result


def find_root_intervals(f, a, b, n=1000):
    intervals, roots = separate_roots(f, a, b, n)
    root_intervals = refine_separation(f, intervals)
    return root_intervals, roots

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

    intervals1, roots1 = find_root_intervals(f, -R, -r)
    intervals2, roots2 = find_root_intervals(f, r, R)
    intervals = intervals1 + intervals2
    roots_cur = roots1 + roots2

    print(f"\n2. Отделение корней методом бисекции:")
    print(f"   Найдено интервалов с корнями: {len(intervals)}")

    for i, interval in enumerate(intervals, 1):
        print(f"   Интервал {i}: [{interval[0]:.6f}, {interval[1]:.6f}]")

    print(f"\n3. Уточнение корней методом бисекции:")
    print("-" * 70)
    print(f"{'№':<3} {'x':<25} {'f(x)':<25} {'Итерации':<10}")
    print("-" * 70)

    roots = []
    for i, interval in enumerate(intervals, 1):
        root, iterations = bisection_method(f, interval[0], interval[1])
        roots.append(root)
        print(f"{i:<3} {root:<25.8f} {f(root):<25.2e} {iterations:<10}")
    print(f"   Найдено точных корней: {len(roots_cur)}")
    for i, root in enumerate(roots_cur, 1):
        print(f"{i:<3} {root:<25.8f} {f(root):<25.2e} {0:<10}")

a3, b3, c3, d3 = 1, 3, 3, 1
analyze_cubic_equation(a3, b3, c3, d3)