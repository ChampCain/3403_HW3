#region imports
import math
#endregion

#region function definitions
def Probability(PDF, args, c, GT=True):
    """
    Calculate the probability that x is greater than or less than c.
    :param PDF: Gaussian PDF function
    :param args: tuple containing (mean, standard deviation)
    :param c: value for which we ask the probability question
    :param GT: boolean deciding if we want probability x > c (True) or x < c (False)
    :return: probability value
    """
    mu, sigma = args
    if GT:
        # P(x > c) = 1 - P(x < c)
        return 1 - Simpson(PDF, (mu - 5 * sigma, c))
    else:
        # P(x < c)
        return Simpson(PDF, (mu - 5 * sigma, c))

def GPDF(args):
    """
    Gaussian Probability Density Function.
    :param args: tuple containing (x, mean, standard deviation)
    :return: value of GPDF at x
    """
    x, mu, sigma = args
    return (1 / (sigma * math.sqrt(2 * math.pi))) * math.exp(-0.5 * ((x - mu) / sigma) ** 2)


def Simpson(fn, limits, N=1000):
    """
    Simpson's 1/3 Rule for numerical integration.
    :param fn: The function to integrate.
    :param limits: A tuple containing (lower_limit, upper_limit).
    :param N: Number of intervals (must be even).
    :return: Integral value.
    """
    a, b = limits  # Unpack lower and upper limits
    if N % 2 != 0:
        N += 1  # Ensure N is even
    h = (b - a) / N
    integral = fn(a) + fn(b)  # Add the endpoints

    for i in range(1, N):
        x = a + i * h
        if i % 2 == 0:
            integral += 2 * fn(x)  # Even intervals
        else:
            integral += 4 * fn(x)  # Odd intervals

    integral *= h / 3
    return integral

def Secant(fcn, x0, x1, maxiter=10, xtol=1e-5):
    """
    Secant Method for finding the root of a function.
    :param fcn: The function for which we want to find the root.
    :param x0: First initial guess.
    :param x1: Second initial guess.
    :param maxiter: Maximum number of iterations.
    :param xtol: Tolerance for convergence.
    :return: A tuple containing the estimated root and the number of iterations performed.
    """
    for i in range(maxiter):
        # Calculate the function values at x0 and x1
        f_x0 = fcn(x0)
        f_x1 = fcn(x1)

        # Avoid division by zero
        if f_x1 - f_x0 == 0:
            break

        # Compute the new estimate using the Secant formula
        x_new = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0)

        # Check for convergence
        if abs(x_new - x1) < xtol:
            return x_new, i + 1  # Return the root and the number of iterations

        # Update x0 and x1 for the next iteration
        x0, x1 = x1, x_new

    # If maxiter is reached, return the current estimate
    return x1, maxiter

def MakeDiagDom(A): # I used AI to create this to make the Guass work
    """
    Make the matrix A diagonally dominant.
    :param A: The input matrix (list of lists).
    :return: A diagonally dominant matrix.
    """
    n = len(A)
    for i in range(n):
        row_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        if abs(A[i][i]) <= row_sum:
            # Find a row to swap with
            for k in range(i + 1, n):
                if abs(A[k][i]) > sum(abs(A[k][j]) for j in range(n) if j != i):
                    A[i], A[k] = A[k], A[i]
                    break
    return A

def GaussSeidel(Aaug, x, Niter=15): # I used AI to create Guass
    """
    Gauss-Seidel method to solve a system of linear equations Ax = b.
    :param Aaug: Augmented matrix [A | b].
    :param x: Initial guess vector.
    :param Niter: Number of iterations.
    :return: The solution vector x.
    """
    Aaug = MakeDiagDom(Aaug)  # Ensure the matrix is diagonally dominant
    n = len(Aaug)
    for _ in range(Niter):
        for i in range(n):
            sigma = sum(Aaug[i][j] * x[j] for j in range(n) if j != i)
            x[i] = (Aaug[i][-1] - sigma) / Aaug[i][i]
    return x

def main():
    '''
    This is a function I created for testing the numerical methods locally.
    :return: None
    '''
    #region testing GPDF
    fx = GPDF((0,0,1))
    print("{:0.5f}".format(fx))  # Does this match the expected value?
    #edregion

    #region testing Simpson
    p=Simpson(GPDF,(0,1,-5,0)) # should return 0.5
    print("p={:0.5f}".format(p))  # Does this match the expected value?
    #endregion

    #region testing Probability
    p1 = Probability(GPDF, (0,1),0,True)
    print("p1={:0.5f}".format(p1))  # Does this match the expected value?
    #endregion
    pass

#endregion

#region function calls
if __name__ == '__main__':
    main()
#endregion