#HW3 part a
#Champ Cain
#MAE 3403

from NumericalMethods import Probability, Secant, GPDF, Simpson

def compute_probability_difference(c, mu, sigma, target_prob, double_sided=False, GT=True):
    """
    Compute the probability difference for Secant method.
    Used to find c that matches target probability.
    """
    if double_sided:
        lower_limit = mu - (c - mu)
        upper_limit = mu + (c - mu)
        prob = Simpson(lambda x: GPDF(x, mu, sigma), (lower_limit, upper_limit))
    else:
        prob = Probability(GPDF, mu, sigma, c, GT=GT)
    return prob - target_prob

def main():
    """
    Main function for interactive user input and probability computation.
    """
    mu = float(input("Enter the mean (μ): "))
    sigma = float(input("Enter the standard deviation (σ): "))
    choice = input(
        "Do you want to specify c and compute P (enter 'c'), or specify P and solve for c (enter 'P')? ").strip().lower()

    if choice == 'c':
        c = float(input("Enter the value of c: "))
        double_sided = input("Do you want a double-sided probability? (y/n): ").strip().lower() == 'y'

        if double_sided:
            lower_limit = mu - (c - mu)
            upper_limit = mu + (c - mu)
            prob = Simpson(lambda x: GPDF(x, mu, sigma), (lower_limit, upper_limit))
            print(f"P({lower_limit:.4f} < x < {upper_limit:.4f} | μ={mu}, σ={sigma}) = {prob:.4f}")
        else:
            GT = input("Do you want P(x > c)? (y/n): ").strip().lower() == 'y'
            prob = Probability(GPDF, mu, sigma, c, GT=GT)
            if GT:
                print(f"P(x > {c:.4f} | μ={mu}, σ={sigma}) = {prob:.4f}")
            else:
                print(f"P(x < {c:.4f} | μ={mu}, σ={sigma}) = {prob:.4f}")

    elif choice == 'p':
        target_prob = float(input("Enter the target probability: "))
        double_sided = input("Do you want a double-sided probability? (y/n): ").strip().lower() == 'y'

        if double_sided:
            def f(c):
                return compute_probability_difference(c, mu, sigma, target_prob, double_sided=True)

            c_initial_guess1 = mu - 2 * sigma
            c_initial_guess2 = mu + 2 * sigma
            c_solution, _ = Secant(f, c_initial_guess1, c_initial_guess2, maxiter=100, xtol=1e-6)

            lower_limit = mu - (c_solution - mu)
            upper_limit = mu + (c_solution - mu)
            print(
                f"To achieve P({lower_limit:.4f} < x < {upper_limit:.4f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")

        else:
            GT = input("Do you want P(x > c)? (y/n): ").strip().lower() == 'y'

            def f(c):
                return compute_probability_difference(c, mu, sigma, target_prob, double_sided=False, GT=GT)

            c_initial_guess1 = mu - 2 * sigma
            c_initial_guess2 = mu + 2 * sigma
            c_solution, _ = Secant(f, c_initial_guess1, c_initial_guess2, maxiter=100, xtol=1e-6)

            if GT:
                print(
                    f"To achieve P(x > {c_solution:.4f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")
            else:
                print(
                    f"To achieve P(x < {c_solution:.4f} | μ={mu}, σ={sigma}) = {target_prob:.4f}, c = {c_solution:.4f}")

    else:
        print("Invalid choice. Please enter 'c' or 'P'.")

if __name__ == "__main__":
    main()