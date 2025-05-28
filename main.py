# Załóżmy, że importujesz wszystkie potrzebne funkcje z pliku leguerre.py
import numpy as np

from leguerre import (
    calculate_c_coeffs,
    calculate_approximating_polynomial_value,  # Twoja "stara" metoda sumy c_k * L_k(x)
    calculate_polynomial_coefficients_from_laguerre,  # Generuje l_kj
    create_coefficients_of_power_polynomial,  # Transformuje c_k -> a_m
    calculate_approximating_polynomial_value_by_horner  # Horner z a_m
)


def main():
    # 1. Zdefiniuj funkcję do aproksymacji f(x)
    def function_to_approximate(x: float) -> float:
        # return x - 1.0  # f(x) = L_1(x)
        # return 1.0 # f(x) = L_0(x)
        # return x**2 - 4*x + 2.0 # f(x) = L_2(x)
         return np.sin(x) + 2

    # 2. Parametry aproksymacji
    n_degree_approximation = 2  # Stopień wielomianu aproksymującego P_N(x)
    num_gauss_nodes_for_c_k = 30  # Liczba węzłów dla kwadratury Gaussa-Laguerre'a

    # Wybierzmy funkcję do opisu w printach
    #func_description = "x - 1.0"
    # func_description = "1.0 (czyli L_0(x))"
    # func_description = "x^2 - 4x + 2.0 (czyli L_2(x))"
    func_description = "sin(x) + 2"

    print(f"\n--- Test Aproksymacji Wielomianami Laguerre'a ---")
    print(f"Funkcja aproksymowana f(x) = {func_description}")
    print(f"Stopień wielomianu aproksymującego N = {n_degree_approximation}")
    print(f"Liczba węzłów Gaussa-Laguerre'a dla c_k = {num_gauss_nodes_for_c_k}")

    # 3. Oblicz współczynniki c_k
    calculated_c_coefficients = calculate_c_coeffs(
        func=function_to_approximate,
        n=n_degree_approximation,
        num_nodes_gauss=num_gauss_nodes_for_c_k
    )

    print("\nObliczone współczynniki c_k [c_0, ..., c_N]:")
    if calculated_c_coefficients:
        for i, c_val in enumerate(calculated_c_coefficients):
            print(f"  c_{i} = {c_val:.8f}")


    print("\nGenerowanie współczynników L_k(x) w bazie potęgowej...")
    all_lk_coeffs = calculate_polynomial_coefficients_from_laguerre(n_degree_approximation)

    # 5. Przekształć c_k na a_m (współczynniki w bazie potęgowej [a_0, ..., a_N])
    print("\nTransformacja współczynników c_k na a_m (dla bazy potęgowej)...")
    a_m_power_coefficients = create_coefficients_of_power_polynomial(
        calculated_c_coefficients,
        all_lk_coeffs
    )
    print(f"\nObliczone współczynniki a_m [a_0, ..., a_N] wielomianu P_N(x) = sum(a_m * x^m):")
    if a_m_power_coefficients:
        for i, a_val in enumerate(a_m_power_coefficients):
            print(f"  a_{i} (wsp. przy x^{i}) = {a_val:.8f}")


    # Wyświetlmy, jak wygląda wielomian P_N(x) w postaci potęgowej
    poly_str_parts = []
    for i, a_val in reversed(list(enumerate(a_m_power_coefficients))):  # Od a_N do a_0
        if abs(a_val) > 1e-9:  # Pokaż tylko niezerowe współczynniki (z tolerancją)
            term = ""
            if a_val > 0 and poly_str_parts:
                term += "+ "
            elif a_val < 0:
                term += "- "

            term += f"{abs(a_val):.4f}"

            if i > 0:
                term += f"x"
                if i > 1:
                    term += f"^{i}"
            poly_str_parts.append(term + " ")

    pn_polynomial_string = "".join(poly_str_parts).strip()
    if not pn_polynomial_string or pn_polynomial_string.startswith("+"):
        pn_polynomial_string = pn_polynomial_string.lstrip("+ ")
    if not pn_polynomial_string:
        pn_polynomial_string = "0.0"

    print(f"\nWielomian aproksymujący P_{n_degree_approximation}(x) ≈ {pn_polynomial_string}")

    # 6. Porównanie wartości P_N(x) z obu metod
    print("\nTestowanie wartości P_N(x) - porównanie metod:")
    test_points_x = [0.0, 0.5, 1.0, 2.0, 5.0]

    print(f"{'x':<8} | {'f(x)':<12} | {'P_N(x) (sum cL_k)':<18} | {'P_N(x) (Horner)':<18} | {'Różnica':<12}")
    print("-" * 75)

    for x_test in test_points_x:
        original_f_value = function_to_approximate(x_test)

        # Metoda 1: Suma c_k * L_k(x)
        approximated_p_n_direct = calculate_approximating_polynomial_value(
            x_eval=x_test,
            c_coefficients=calculated_c_coefficients
        )

        # Metoda 2: Horner z a_m
        approximated_p_n_horner = calculate_approximating_polynomial_value_by_horner(
            x=x_test,
            a_coefficients=a_m_power_coefficients
        )

        diff_methods = abs(approximated_p_n_direct - approximated_p_n_horner)

        print(
            f"{x_test:<8.2f} | {original_f_value:<12.6f} | {approximated_p_n_direct:<18.6f} | {approximated_p_n_horner:<18.6f} | {diff_methods:<12.2e}")

    print("\n--- Koniec Testu ---")


if __name__ == '__main__':
    main()