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

a = 0.6
b = 1.4
N = 20 # Для Симпсона четное

trap_result = composite_trapezoidal(f, a, b, N)
simp_result = composite_simpson(f, a, b, N)

print(f"Составная формула трапеций: {trap_result}")
print(f"Составная формула Симпсона: {simp_result}")

print(f"С помощью калькулятора: 0.2222658545486588")