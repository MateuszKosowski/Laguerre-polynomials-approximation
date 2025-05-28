import numpy as np
import math
from leguerre import calculate_approximating_polynomial_value, calculate_c_coeffs

def main():
    # 1. Zdefiniuj funkcję do aproksymacji f(x)
    def function_to_approximate(x: float) -> float:
        return x - 1.0
        # return x**2 - 6*x + 5.0
        # return np.exp(-x/2.0)

    # 2. Parametry aproksymacji
    n_degree_approximation = 3  # Stopień wielomianu aproksymującego P_N(x)
    num_gauss_nodes_for_c_k = 20  # Liczba węzłów dla kwadratury Gaussa-Laguerre'a przy liczeniu c_k

    print(f"\n--- Test Aproksymacji Wielomianami Laguerre'a ---")
    print(f"Funkcja aproksymowana f(x) = x - 1")
    print(f"Stopień wielomianu aproksymującego N = {n_degree_approximation}")
    print(f"Liczba węzłów Gaussa-Laguerre'a dla c_k = {num_gauss_nodes_for_c_k}")

    # 3. Współczynniki c_k
    calculated_c_coefficients = calculate_c_coeffs(
        func=function_to_approximate,
        n=n_degree_approximation,
        num_nodes_gauss=num_gauss_nodes_for_c_k
    )

    print("\nObliczone współczynniki c_k:")
    for i, c_val in enumerate(calculated_c_coefficients):
        print(f"  c_{i} = {c_val:.8f}")


    # 4. Wartości wielomianu aproksymującego P_N(x) w kilku punktach

    print("\nTestowanie wartości wielomianu aproksymującego P_N(x):")
    test_points_x = [0.0, 0.5, 1.0, 2.0, 5.0]

    for x_test in test_points_x:
        original_f_value = function_to_approximate(x_test)
        approximated_P_N_value = calculate_approximating_polynomial_value(
            x_eval=x_test,
            c_coefficients=calculated_c_coefficients
        )

        error = abs(original_f_value - approximated_P_N_value)
        print(f"  Dla x = {x_test:.2f}: f(x) = {original_f_value:.6f}, P_N(x) = {approximated_P_N_value:.6f}, Błąd = {error:.2e}")

if __name__ == '__main__':
    main()
