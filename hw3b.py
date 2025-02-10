#HW3 part b
#Champ Cain
#MAE 3403

import math
from NumericalMethods import Simpson  # Import the updated Simpson function

def gamma_function(alpha):
    """
    Compute the gamma function using math.gamma.
    :param alpha: Argument to the gamma function.
    :return: Value of the gamma function.
    """
    return math.gamma(alpha)

def compute_Km(m):
    """
    Compute the normalization constant K_m for the t-distribution.
    :param m: Degrees of freedom.
    :return: Value of K_m.
    """
    numerator = gamma_function((m + 1) / 2)
    denominator = math.sqrt(m * math.pi) * gamma_function(m / 2)
    return numerator / denominator

def t_distribution_pdf(u, m):
    """
    Compute the PDF of the t-distribution at a given point u.
    :param u: Point at which to evaluate the PDF.
    :param m: Degrees of freedom.
    :return: Value of the PDF at u.
    """
    return (1 + (u**2) / m) ** (-(m + 1) / 2)

def t_distribution_cdf(z, m, N=1000):
    """
    Compute the CDF of the t-distribution using numerical integration.
    :param z: Upper limit of integration.
    :param m: Degrees of freedom.
    :param N: Number of intervals for numerical integration.
    :return: Value of the CDF at z.
    """
    Km = compute_Km(m)
    # Integrate from -infinity to z. For practical purposes, integrate from -large_value to z.
    lower_limit = -100  # Approximating -infinity
    upper_limit = z
    # Define the integrand as a function of u
    integrand = lambda u: t_distribution_pdf(u, m)
    # Compute the integral using Simpson's rule
    integral = Simpson(integrand, (lower_limit, upper_limit), N=N)
    return Km * integral

def main():
    """
    Main function to interact with the user and compute the t-distribution CDF.
    """
    # Get user input
    m = int(input("Enter the degrees of freedom (m): "))
    z = float(input("Enter the value of z: "))

    # Compute the CDF
    cdf_value = t_distribution_cdf(z, m)

    # Print the result
    print(f"F(z) = {cdf_value:.6f} for m = {m} and z = {z}")

if __name__ == "__main__":
    main()