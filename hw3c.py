#HW3 part c
#Champ Cain
#MAE 3403

#Region imports
import math



#I had to use some AI help with the "IF" and "ELSE" statements in setting them up correctly

def is_symmetric(A):
    """
    Check if a matrix is symmetric.
    :param A: The matrix to check.
    :return: True if symmetric, False otherwise.
    """
    n = len(A)
    for i in range(n):
        for j in range(n):
            if A[i][j] != A[j][i]:
                return False
    return True

def is_positive_definite(A):
    """
    Check if a matrix is positive definite using the Cholesky decomposition method.
    :param A: The matrix to check.
    :return: True if positive definite, False otherwise.
    """
    try:
        cholesky_decomposition(A)
        return True
    except ValueError:
        return False

def cholesky_decomposition(A):
    """
    Perform Cholesky decomposition on a symmetric, positive definite matrix.
    :param A: The matrix to decompose.
    :return: The lower triangular matrix L.
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1):
            sum_ = sum(L[i][k] * L[j][k] for k in range(j))
            if i == j:
                if A[i][i] - sum_ <= 0:
                    raise ValueError("Matrix is not positive definite")
                L[i][j] = math.sqrt(A[i][i] - sum_)
            else:
                L[i][j] = (A[i][j] - sum_) / L[j][j]
    return L

def doolittle_decomposition(A):
    """
    Perform Doolittle LU decomposition on a matrix.
    :param A: The matrix to decompose.
    :return: The lower triangular matrix L and the upper triangular matrix U.
    """
    n = len(A)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        # Upper triangular matrix
        for k in range(i, n):
            sum_ = sum(L[i][j] * U[j][k] for j in range(i))
            U[i][k] = A[i][k] - sum_

        # Lower triangular matrix
        for k in range(i, n):
            if i == k:
                L[i][i] = 1.0  # Diagonal elements are 1
            else:
                sum_ = sum(L[k][j] * U[j][i] for j in range(i))
                L[k][i] = (A[k][i] - sum_) / U[i][i]

    return L, U

def solve_lower_triangular(L, b):
    """
    Solve the system Ly = b where L is a lower triangular matrix.
    :param L: The lower triangular matrix.
    :param b: The right-hand side vector.
    :return: The solution vector y.
    """
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        y[i] = b[i] - sum(L[i][j] * y[j] for j in range(i))
        y[i] /= L[i][i]
    return y

def solve_upper_triangular(U, y):
    """
    Solve the system Ux = y where U is an upper triangular matrix.
    :param U: The upper triangular matrix.
    :param y: The right-hand side vector.
    :return: The solution vector x.
    """
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] /= U[i][i]
    return x

def solve_cholesky(A, b):
    """
    Solve the system Ax = b using Cholesky decomposition.
    :param A: The coefficient matrix.
    :param b: The right-hand side vector.
    :return: The solution vector x.
    """
    L = cholesky_decomposition(A)
    y = solve_lower_triangular(L, b)
    x = solve_upper_triangular([list(row) for row in zip(*L)], y)  # Transpose L to get L.T
    return x

def solve_doolittle(A, b):
    """
    Solve the system Ax = b using Doolittle LU decomposition.
    :param A: The coefficient matrix.
    :param b: The right-hand side vector.
    :return: The solution vector x.
    """
    L, U = doolittle_decomposition(A)
    y = solve_lower_triangular(L, b)
    x = solve_upper_triangular(U, y)
    return x

def main():
    """
    Main function to solve the given systems of equations.
    """
    # Problem 1
    A1 = [
        [1, -1, 3, 2],
        [-1, 5, -5, -2],
        [3, -5, 19, 3],
        [2, -2, 3, 21]
    ]
    b1 = [15, -35, 94, 1]

    if is_symmetric(A1) and is_positive_definite(A1):
        x1 = solve_cholesky(A1, b1)
        method1 = "Cholesky"
    else:
        x1 = solve_doolittle(A1, b1)
        method1 = "Doolittle"

    print(f"Solution to Problem 1 (using {method1} method):")
    print(x1)

    # Problem 2
    A2 = [
        [4, 2, 4, 0],
        [2, 2, 3, 2],
        [4, 3, 6, 3],
        [0, 2, 3, 9]
    ]
    b2 = [20, 36, 60, 122]

    if is_symmetric(A2) and is_positive_definite(A2):
        x2 = solve_cholesky(A2, b2)
        method2 = "Cholesky"
    else:
        x2 = solve_doolittle(A2, b2)
        method2 = "Doolittle"

    print(f"Solution to Problem 2 (using {method2} method):")
    print(x2)

if __name__ == "__main__":
    main()