import numpy

def f(x):
    return numpy.cos(x) / (x + 1)

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

def runge_rule(J_h, J_2h, p):
    Chp = (J_h - J_2h) / (2**p - 1)
    J = J_h + Chp

    return Chp, J

a = 0.6
b = 1.4
N = 20 # Для Симпсона четное

trap_result = composite_trapezoidal(f, a, b, N)
simp_result = composite_simpson(f, a, b, N)
trap_result_2h = composite_trapezoidal(f, a, b, N // 2)
simp_result_2h = composite_simpson(f, a, b, N // 2)

trap_error, trap_exact = runge_rule(trap_result, trap_result_2h, 2)
simp_error, simp_exact = runge_rule(simp_result, simp_result_2h, 4)

print(f"Формула трапеций:")
print(f"    J_h: {trap_result}")
print(f"    J_2h: {trap_result_2h}")
print(f"    Точное решение: {trap_exact}\n")
print(f"    Погрешность по Рунге: {trap_error}")
print(f"    Погрешность J - J_h: {trap_exact - trap_result}\n")


print(f"Формула Симпсона:")
print(f"    J_h: {simp_result}")
print(f"    J_2h: {simp_result_2h}")
print(f"    Точное решение: {simp_exact}\n")
print(f"    Погрешность по Рунге: {simp_error}")
print(f"    Погрешность J - J_h: {simp_exact - simp_result}\n")
