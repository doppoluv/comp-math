import math

a = 0.6
b = 1.4
N = 20

def f(x):
    return math.cos(x) / (x + 1)

def composite_trapezoidal(f, a, b, N):
    h = (b - a) / N
    x = [a + i * h for i in range(N + 1)]
    fx = [f(xi) for xi in x]
    
    integral = h / 2 * (fx[0] + 2 * sum(fx[1:N]) + fx[N])
    return integral

def composite_simpson(f, a, b, N):
    h = (b - a) / N
    x = [a + i * h for i in range(N + 1)]
    fx = [f(xi) for xi in x]
    
    integral = h / 3 * (fx[0] + 4 * sum(fx[1:N:2]) + 2 * sum(fx[2:N:2]) + fx[N])
    return integral

def find_sign_changes(f, a, b, N):
    h = (b - a) / N
    x = [a + i * h for i in range(N + 1)]
    fx = [f(xi) for xi in x]
    segments = []
    start = a
    start_sign = 1 if fx[0] >= 0 else -1

    for i in range(1, len(fx)):
        current_sign = 1 if fx[i] >= 0 else -1
        if current_sign != start_sign or i == len(fx) - 1:
            end = x[i] if i == len(fx) - 1 else x[i - 1]
            segments.append((start, end))
            start = x[i - 1]
            start_sign = current_sign
    return segments

def integrate_with_sign_check(f, a, b, N, method):
    segments = find_sign_changes(f, a, b, N)
    integral = 0
    for seg_a, seg_b in segments:
        seg_N = max(2, int(N * (seg_b - seg_a) / (b - a)))
        if method == 'trap':
            integral += composite_trapezoidal(f, seg_a, seg_b, seg_N)
        elif method == 'simp':
            integral += composite_simpson(f, seg_a, seg_b, seg_N)
    return integral

def runge_rule(J1, J2, k, p):
    Chp = (J1 - J2) / (1 - 1 / k**p)
    J = J1 + Chp
    return Chp, J

def p_counting(J1, J2, J3, k):
    try:
        return math.log((abs(J1 - J2) / abs(J2 - J3)), k)
    except (ValueError, ZeroDivisionError):
        return float('nan')



trap_result = integrate_with_sign_check(f, a, b, N, 'trap')
simp_result = integrate_with_sign_check(f, a, b, N, 'simp')
trap_result_h_2 = integrate_with_sign_check(f, a, b, 2 * N, 'trap')
simp_result_h_2 = integrate_with_sign_check(f, a, b, 2 * N, 'simp')
trap_result_h_4 = integrate_with_sign_check(f, a, b, 4 * N, 'trap')
simp_result_h_4 = integrate_with_sign_check(f, a, b, 4 * N, 'simp')

trap_error, trap_exact = runge_rule(trap_result, trap_result_h_2, 2, 2)
simp_error, simp_exact = runge_rule(simp_result, simp_result_h_2, 2, 4)

p_trap = p_counting(trap_result, trap_result_h_2, trap_result_h_4, 2)
P_simp = p_counting(simp_result, simp_result_h_2, simp_result_h_4, 2)
trap_error_2, trap_exact_2 = runge_rule(trap_result, trap_result_h_2, 2, p_trap)
simp_error_2, simp_exact_2 = runge_rule(simp_result, simp_result_h_2, 2, P_simp)

print(f"Формула трапеций:")
print(f"    J_h: {trap_result}")
print(f"  p = 2:")
print(f"    Точное решение: {trap_exact}")
print(f"    Погрешность: {trap_error}")
print(f"  p = {p_trap}:")
print(f"    Точное решение: {trap_exact_2}")
print(f"    Погрешность: {trap_error_2}\n")

print(f"Формула Симпсона:")
print(f"    J_h: {simp_result}")
print(f"  p = 4:")
print(f"    Точное решение: {simp_exact}")
print(f"    Погрешность: {simp_error}")
print(f"  p = {P_simp}:")
print(f"    Точное решение: {simp_exact_2}")
print(f"    Погрешность: {simp_error_2}\n")

segments = find_sign_changes(f, a, b, N)
print("Подотрезки:", [(seg_a, seg_b) for seg_a, seg_b in segments])