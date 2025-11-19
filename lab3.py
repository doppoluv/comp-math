import numpy as np
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

def lu_decomposition(A):
    n = A.shape[0]
    A = A.copy()
    L = np.eye(n)

    for k in range(n-1):
        if abs(A[k, k]) < 1e-12:
            raise ValueError(f"A[{k},{k}] = {A[k,k]}")

        for i in range(k+1, n):
            l_ik = A[i, k] / A[k, k]
            L[i, k] = l_ik
            
            print(f"l_{i+1}{k+1} = a_{i+1}{k+1} / a_{k+1}{k+1} = {A[i,k]:.3f} / {A[k,k]:.3f} = {l_ik:.3f}")
            
            for j in range(k, n):
                A[i, j] = A[i, j] - l_ik * A[k, j]

        print(f"  Матрица после шага {k+1}:")
        print(np.array2string(A, precision=3, suppress_small=True))
        print()

    U = A.copy()
    return L, U


# Ly = b
def forward_substitution(L, b):
    n = len(b)
    y = np.zeros(n)

    for i in range(n):
        y[i] = (b[i] - sum(L[i, j] * y[j] for j in range(i))) / L[i, i]

    return y


# Ux = y
def backward_substitution(U, y):
    n = len(y)
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i, j] * x[j] for j in range(i + 1, n))) / U[i, i]

    return x


def jacobi(A, b, eps=1e-8, max_iter=10000):
    n = len(b)
    x = np.zeros(n)
    
    for k in range(max_iter):
        x_new = np.zeros(n)
        
        for i in range(n):
            sum_j = sum(A[i, j] * x[j] for j in range(n) if j != i)
            x_new[i] = (b[i] - sum_j) / A[i, i]
        
        if np.linalg.norm(x_new - x) < eps:
            return x_new
        
        x = x_new
    
    print("Максимальное число итераций в Якоби")
    return x


def gauss_seidel(A, b, eps=1e-8, max_iter=10000):
    n = len(b)
    x = np.zeros(n)
    
    for k in range(max_iter):
        x_old = x.copy()
        
        for i in range(n):
            sum1 = sum(A[i, j] * x[j] for j in range(i))
            sum2 = sum(A[i, j] * x_old[j] for j in range(i + 1, n))
            x[i] = (b[i] - sum1 - sum2) / A[i, i]

        if np.linalg.norm(x - x_old) < eps:
            return x
    
    print("Максимальное число итераций в Гауссе-Зейделе")
    return x


def thomas_algorithm(A, b):
    a = np.diag(A, -1).copy()
    d = np.diag(A).copy()
    c = np.diag(A, 1).copy()
    n = len(b)

    # Прямой ход
    alpha = np.zeros(n - 1)
    beta = np.zeros(n)
    alpha[0] = -c[0] / d[0]
    beta[0] = b[0] / d[0]
    for i in range(1, n):
        denom = d[i] + a[i - 1] * alpha[i - 1]
        if i < n - 1:
            alpha[i] = -c[i] / denom
        beta[i] = (b[i] - a[i - 1] * beta[i - 1]) / denom

    # Обратный ход
    x = np.zeros(n)
    x[-1] = beta[-1]
    for i in range(n - 2, -1, -1):
        x[i] = beta[i] + alpha[i] * x[i + 1]

    return x


def check_solution(x, A, b):
    residual = np.dot(A, x) - b
    return np.linalg.norm(residual)



A = np.array([
    [8.2, -3.2, 14.2, 14.8],
    [5.6, -12.0, 15.0, -6.4],
    [5.7, 3.6, -12.4, -2.3],
    [6.8, 13.2, -6.3, -8.7]
], dtype=float)
b = np.array([-8.4, 4.5, 3.3, 14.3], dtype=float)
n = A.shape[0]

L, U = lu_decomposition(A)
y = forward_substitution(L, b)
x_lu = backward_substitution(U, y)
x_jacobi = jacobi(A, b)
x_gs = gauss_seidel(A, b)
# x_thomas = thomas_algorithm(A, b)

print("\nLU разложение")
print("   Матрица L:\n", np.array2string(L, precision=3, suppress_small=True))
print("   Матрица U:\n", np.array2string(U, precision=3, suppress_small=True))

print("\nНормы и числа обусловленности:")
norm_types = [1, 2, np.inf]
for k in norm_types:
    norm_A = np.linalg.norm(A, k)
    cond_A = np.linalg.cond(A, k)
    norm_L = np.linalg.norm(L, k)
    cond_L = np.linalg.cond(L, k)
    norm_U = np.linalg.norm(U, k)
    cond_U = np.linalg.cond(U, k)

    print(f"   {k}:")
    print(f"     ||A||_{k} = {norm_A:.6f},  ν_{k}(A) = {cond_A:.6f}")
    print(f"     ||L||_{k} = {norm_L:.6f},  ν_{k}(L) = {cond_L:.6f}")
    print(f"     ||U||_{k} = {norm_U:.6f},  ν_{k}(U) = {cond_U:.6f}")

print("\nРешения системы Ax = b:")
print(f"   LU-разложение:  {x_lu}")
print(f"   Якоби:          {x_jacobi}")
print(f"   Гаусс-Зейдель:  {x_gs}")
# print(f"   Прогонка:       {x_thomas}")

print("\nПроверка решений (||Ax - b||):")
print(f"   LU:             {check_solution(x_lu, A, b):.2e}")
print(f"   Якоби:          {check_solution(x_jacobi, A, b):.2e}")
print(f"   Гаусс-Зейдель:  {check_solution(x_gs, A, b):.2e}")
# print(f"   Прогонка:     {check_solution(x_thomas, A, b):.2e}")

x_exact = np.linalg.solve(A, b)
print(f"\nТочное решение: {x_exact}")