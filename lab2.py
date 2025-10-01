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

def runge_rule(J1, J2, k, p):
    Chp = (J1 - J2) / (1 - 1 / k**p)
    J = J1 + Chp

    return Chp, J

def p_counting(J1, J2, J3, k):
    return math.log(((J1 - J2) / (J2 - J3)), k)



trap_result = composite_trapezoidal(f, a, b, N)
simp_result = composite_simpson(f, a, b, N)
trap_result_h_2 = composite_trapezoidal(f, a, b, 2 * N)
simp_result_h_2 = composite_simpson(f, a, b, 2 * N)
trap_result_h_4 = composite_trapezoidal(f, a, b, 4 * N)
simp_result_h_4 = composite_simpson(f, a, b, 4 * N)

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
