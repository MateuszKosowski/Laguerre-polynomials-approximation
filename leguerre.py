import numpy as np

# Aproksymacja to jest proces, w którym staramy się znaleźć funkcję (wielomian), która najlepiej pasuje do danej funkcji na określonym przedziale.

# Wzór wielomianu aproksymacyjnego Laguerre'a:
# L_n(x) = sum(c_k * L_k(x)), gdzie c_k to współczynniki aproksymacji

def calculate_laguerre_polynomial_value(n: int, x: float) -> float:
    # n - stopień wielomianu
    # x - punkt, w którym obliczamy wielomian
    # Zwracamy wartość wielomianu Laguerre'a L_k(x)  dla danego n i x


    if n == 0:
        return 1.0
    elif n == 1:
        return x - 1
    else:
        leguerre_prev_prev = 1.0  # Wielomian stopnia 0
        leguerre_prev = x - 1.0  #  Wielomian stopnia 1

        for i in range(2, n + 1):

            # Iteracyjne obliczanie wielomianu Laguerre'a n stopnia za pomocą wzoru rekurencyjnego
            # Wzór rekurencyjny z wykładu:
            leguerre_current = ((x - 2*(i-1) -1) * leguerre_prev - ((i-1)*(i-1)) * leguerre_prev_prev)

            # Aktualizujemy poprzednie wielomiany
            leguerre_prev_prev, leguerre_prev = leguerre_prev, leguerre_current

        return leguerre_prev


# Zakładamy, że funkcja func jest zdefiniowana na przedziale [0, +niekończoność) i, że jest ortogonalna dla wielomianów Laguerre'a.
# Ortogonalność funkcji oznacza, że iloczyn skalarny dwóch różnych funkcji jest równy zeru.
# Czyli: całka na przedziale [0, +∞) z iloczynu dwóch różnych funkcji jest równa zeru.

# c_k = ∫_0^∞ f(x) * L_k(x) * e^(-x) dx / ∫_0^∞ e^(-x) * L_k(x)^2 dx

def calculate_c_coeffs(
        func: callable,
        n: int,
        num_nodes_gauss: int = 20,
    ) -> list[float]:

    c_k_list =  [0.0] * (n + 1) # n + 1, ponieważ mamy współczynniki od 0 do n

    # Obliczamy węzły i ich wagi kwadratury Gaussa-Laguerre'a
    x_gauss_nodes, w_gauss_weights = np.polynomial.laguerre.laggauss(num_nodes_gauss)

    for k in range(n + 1):

        # Obliczamy tą całkę na górze wzoru za pomocą kwadratury Gaussa-Laguerre'a
        integral_numerator_val = 0.0
        for i_node in range(num_nodes_gauss):
            current_node_x = x_gauss_nodes[i_node]
            current_weight = w_gauss_weights[i_node]

            val_f_at_node = func(current_node_x)
            val_leguerre_k_at_node = calculate_laguerre_polynomial_value(k, current_node_x)

            integral_numerator_val += current_weight * val_f_at_node * val_leguerre_k_at_node

        # Obliczamy tą całkę na dole wzoru za pomocą kwadratury Gaussa-Laguerre'a
        integral_denominator_val = 0.0
        for i_node in range(num_nodes_gauss):
            current_node_x = x_gauss_nodes[i_node]
            current_weight = w_gauss_weights[i_node]

            val_leguerre_k_at_node = calculate_laguerre_polynomial_value(k, current_node_x)

            integral_denominator_val += current_weight * (val_leguerre_k_at_node ** 2)

        c_k_list[k] = integral_numerator_val / integral_denominator_val

    return c_k_list

# Funkcja ta oblicza wartość wielomianu aproksymacyjnego w punkcie x
def calculate_approximating_polynomial_value(
        x_eval: float,
        c_coefficients: list[float]
    ) -> float:

    n = len(c_coefficients) - 1

    pn_x_value = 0.0
    for k in range(n + 1):
        c_k = c_coefficients[k]
        l_k_at_x_eval = calculate_laguerre_polynomial_value(k, x_eval)
        pn_x_value += c_k * l_k_at_x_eval

    return pn_x_value

