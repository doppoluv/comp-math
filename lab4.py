import numpy as np

def df(f, x, h=1e-8):
    return (f(x + h) - f(x - h)) / (2 * h)

def separate_roots(f, a, b, n=100):
    x = np.linspace(a, b, n + 1)
    intervals = []
    
    # Ищем отрезки со сменой знака
    for i in range(n):
        if f(x[i]) * f(x[i + 1]) < 0:
            intervals.append([x[i], x[i + 1]])
    
    return intervals

def refine_separation(f, intervals, max_iter=100):
    result = []
    
    for interval in intervals:
        a, b = interval
        x_test = np.linspace(a, b, 20)
        df_vals = [df(f, x) for x in x_test]
        
        # Если производная сохраняет знак - корень единственный
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

def find_root_intervals(f, a=0, b=1, n=100):
    intervals = separate_roots(f, a, b, n)
    root_intervals = refine_separation(f, intervals)
    return root_intervals

def f(x):
    return np.cos(x)

a = -10
b = 10

root_intervals = find_root_intervals(f, a, b, 100)
print(f"\nНайдено отрезков: {len(root_intervals)}\n")
for i, interval in enumerate(root_intervals, 1):
    a_i, b_i = interval
    print(f"Корень {i}:")
    print(f"  Отрезок: [{a_i}, {b_i}]")
    print(f"  f(a) = {f(a_i)}, f(b) = {f(b_i)}\n")
