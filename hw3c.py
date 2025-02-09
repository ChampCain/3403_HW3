import numpy as np

#I had to use some AI help with the "IF" and "ELSE" statements in setting them up correctly

def is_symmetric(A):
    """
    Check if a matrix is symmetric.
    :param A: The matrix to check.
    :return: True if symmetric, False otherwise.
    """
    return np.allclose(A, A.T)

def is_positive_definite(A):
    """
    Check if a matrix is positive definite.
    :param A: The matrix to check.
    :return: True if positive definite, False otherwise.
    """
    try:
        np.linalg.cholesky(A)
        return True
    except np.linalg.LinAlgError:
        return False

def cholesky_decomposition(A):
    """
    Perform Cholesky decomposition on a symmetric, positive definite matrix.
    :param A: The matrix to decompose.
    :return: The lower triangular matrix L.
    """
    return np.linalg.cholesky(A)

def doolittle_decomposition(A):
    """
    Perform Doolittle LU decomposition on a matrix.
    :param A: The matrix to decompose.
    :return: The lower triangular matrix L and the upper triangular matrix U.
    """
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))

    for i in range(n):
        # Upper triangular matrix
        for k in range(i, n):
            sum_ = sum(L[i][j] * U[j][k] for j in range(i))
            U[i][k] = A[i][k] - sum_

        # Lower triangular matrix
        for k in range(i, n):
            if i == k:
                L[i][i] = 1  # Diagonal elements are 1
            else:
                sum_ = sum(L[k][j] * U[j][i] for j in range(i))
                L[k][i] = (A[k][i] - sum_) / U[i][i]

    return L, U

def solve_cholesky(A, b):
    """
    Solve the system Ax = b using Cholesky decomposition.
    :param A: The coefficient matrix.
    :param b: The right-hand side vector.
    :return: The solution vector x.
    """
    L = cholesky_decomposition(A)
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(L.T, y)
    return x

def solve_doolittle(A, b):
    """
    Solve the system Ax = b using Doolittle LU decomposition.
    :param A: The coefficient matrix.
    :param b: The right-hand side vector.
    :return: The solution vector x.
    """
    L, U = doolittle_decomposition(A)
    y = np.linalg.solve(L, b)
    x = np.linalg.solve(U, y)
    return x

def main():
    """
    Main function to solve the given systems of equations.
    """
    # Problem 1
    A1 = np.array([
        [1, -1, 3, 2],
        [-1, 5, -5, -2],
        [3, -5, 19, 3],
        [2, -2, 3, 21]
    ])
    b1 = np.array([15, -35, 94, 1])

    if is_symmetric(A1) and is_positive_definite(A1):
        x1 = solve_cholesky(A1, b1)
        method1 = "Cholesky"
    else:
        x1 = solve_doolittle(A1, b1)
        method1 = "Doolittle"

    print(f"Solution to Problem 1 (using {method1} method):")
    print(x1)

    # Problem 2
    A2 = np.array([
        [4, 2, 4, 0],
        [2, 2, 3, 2],
        [4, 3, 6, 3],
        [0, 2, 3, 9]
    ])
    b2 = np.array([20, 36, 60, 122])

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