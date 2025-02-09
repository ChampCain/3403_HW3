from NumericalMethods import Probability, Secant, GPDF

def compute_probability_difference(c, mu, sigma, target_prob, double_sided=False, GT=True):
    """
    Helper function to compute the difference between the target probability and the computed probability.
    This function is used by the Secant method to find the value of c that matches the target probability.
    :param c: The value of c.
    :param mu: The mean of the normal distribution.
    :param sigma: The standard deviation of the normal distribution.
    :param target_prob: The target probability.
    :param double_sided: Whether to compute a double-sided probability.
    :param GT: Whether to compute P(x > c) or P(x < c).
    :return: The difference between the target probability and the computed probability.

    I had some help with AI on this one. Below are the instructions I put together to iron out the complexities.

    Enter the mean (μ): 100
    Enter the standard deviation (σ): 12.5
    Do you want to specify c and compute P (enter 'c'), or specify P and solve for c (enter 'P')? P
    Enter the target probability: 0.95
    Do you want a double-sided probability? (y/n): n
    Do you want P(x > c)? (y/n): n


    """
    if double_sided:
        # Double-sided probability: P(mu - (c - mu) < x < mu + (c - mu))
        lower_limit = mu - (c - mu)
        upper_limit = mu + (c - mu)
        prob = Probability(GPDF, (mu, sigma), upper_limit, GT=False) - Probability(GPDF, (mu, sigma), lower_limit, GT=False)
    else:
        # Single-sided probability: P(x < c) or P(x > c)
        prob = Probability(GPDF, (mu, sigma), c, GT=GT)
    return prob - target_prob

def main():
    """
    Main function to interact with the user and compute probabilities or solve for c.
    """
    # Get user input
    mu = float(input("Enter the mean (μ): "))
    sigma = float(input("Enter the standard deviation (σ): "))
    choice = input("Do you want to specify c and compute P (enter 'c'), or specify P and solve for c (enter 'P')? ").strip().lower()

    if choice == 'c':
        # User specifies c and wants to compute P
        c = float(input("Enter the value of c: "))
        double_sided = input("Do you want a double-sided probability? (y/n): ").strip().lower() == 'y'
        if double_sided:
            # Double-sided probability
            lower_limit = mu - (c - mu)
            upper_limit = mu + (c - mu)
            prob = Probability(GPDF, (mu, sigma), upper_limit, GT=False) - Probability(GPDF, (mu, sigma), lower_limit, GT=False)
            print(f"P({lower_limit:.2f} < x < {upper_limit:.2f} | μ={mu}, σ={sigma}) = {prob:.4f}")
        else:
            # Single-sided probability
            GT = input("Do you want P(x > c)? (y/n): ").strip().lower() == 'y'
            prob = Probability(GPDF, (mu, sigma), c, GT=GT)
            if GT:
                print(f"P(x > {c:.2f} | μ={mu}, σ={sigma}) = {prob:.4f}")
            else:
                print(f"P(x < {c:.2f} | μ={mu}, σ={sigma}) = {prob:.4f}")

    elif choice == 'p':
        # User specifies P and wants to solve for c
        target_prob = float(input("Enter the target probability: "))
        double_sided = input("Do you want a double-sided probability? (y/n): ").strip().lower() == 'y'
        if double_sided:
            # Double-sided probability: solve for c such that P(mu - (c - mu) < x < mu + (c - mu)) = target_prob
            # Use the Secant method to find c
            def f(c):
                return compute_probability_difference(c, mu, sigma, target_prob, double_sided=True)
            c_initial_guess = mu + sigma  # Initial guess for c
            c_solution, _ = Secant(f, c_initial_guess, c_initial_guess + 1, maxiter=100, xtol=1e-6)
            lower_limit = mu - (c_solution - mu)
            upper_limit = mu + (c_solution - mu)
            print(f"To achieve P({lower_limit:.2f} < x < {upper_limit:.2f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")
        else:
            # Single-sided probability: solve for c such that P(x < c) or P(x > c) = target_prob
            GT = input("Do you want P(x > c)? (y/n): ").strip().lower() == 'y'
            def f(c):
                return compute_probability_difference(c, mu, sigma, target_prob, double_sided=False, GT=GT)
            c_initial_guess = mu + sigma  # Initial guess for c
            c_solution, _ = Secant(f, c_initial_guess, c_initial_guess + 1, maxiter=100, xtol=1e-6)
            if GT:
                print(f"To achieve P(x > {c_solution:.4f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")
            else:
                print(f"To achieve P(x < {c_solution:.4f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")

    else:
        print("Invalid choice. Please enter 'c' or 'P'.")

if __name__ == "__main__":
    main()