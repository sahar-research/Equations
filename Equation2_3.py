import sympy as sp
import time

start_time = time.time()

N = 2
K = 5
L = 3.6
M = 60
b = 0.6

#definition parameter
a_ij = sp.Matrix([[sp.symbols(f'a_{i}_{2*m-2}') for m in range(1, K + 1)] for i in range(1, 2 * N)])

# تنظیم مقادیر a_{i,0} برابر صفر
for i in range(1, 2 * N):
    a_ij[i - 1, 0] = 0


def compute_C_n_2m_2(n, m, M, L, N, b, i):
    C_n_2m_2 = sp.Sum((2 / (M * L)) * a_ij[sp.symbols('i')-1, m-1] * sp.sin(n * sp.symbols('i') * b * sp.pi / L),
                      (sp.symbols('i'), 1, 2 * N - 1)).doit()
    return C_n_2m_2

equations = []


for k in range(1, K+1):

    for l in range(1, 2 * N):
        result_l = 0
        for n in range(1, 4):

            omega_n = 24.62 * n ** 2 * 3.14**2
            for m in range(1, k+1):  # m به عنوان شمارنده

                C_n_2m_2 = compute_C_n_2m_2(n, m, M, L, N, b, l)

                term = (-1)**(k + m) * C_n_2m_2 * omega_n ** (2*(k - m)) * sp.factorial(2 * m - 2)
                term *= sp.sin((n * sp.pi * l) / (2 * N))
                result_l += term
        equations.append(sp.Eq(result_l, 0))

valid_equations = [eq for eq in equations if isinstance(eq, sp.Equality)]

rounded_equations = [sp.Eq(sp.N(eq.lhs, 2), sp.N(eq.rhs, 2)) for eq in valid_equations]

remaining_vars = a_ij[:, 1:]
solution = sp.solve(rounded_equations, remaining_vars)

sp.pprint(solution)

with open('aij_solutions.txt', 'w') as file:
    file.write("Solutions for a_ij (with coefficients rounded to 2 decimals):\n")
    for sol in solution:
        file.write(f'{sol} = {solution[sol]}\n')

end_time = time.time()
execution_time = end_time - start_time

with open('aij_solutions.txt', 'a') as file:  # 'a' برای الحاق به فایل
    file.write("\nExecution Time: {:.4f} seconds\n".format(execution_time))



