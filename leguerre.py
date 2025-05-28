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

# TODO: Przekształcić P_N(x) = Σ (c_k*L_k(x)) do postaci potęgowej P_N(x) = Σ a_m * x^m.
# Każdy wielomian Laguerre'a L_k(x) można wyrazić jako wielomian w postaci potęgowej:
# L_k(x) = l_{k,0} x^0 + l_{k,1} x^1 + l_{k,2} x^2 + ... + l_{k,k} x^k
    # Aby znaleźć współczynniki wielomianu L_i(x) na podstawie współczynników L_{i-1}(x) i L_{i-2}(x), wykonujemy następujące operacje na wielomianach:
    # 1) Mnożenie L_{i-1}(x) przez x
    # 2) Mnożenie L_{i-1}(x) przez skalar -(2i-1)
    # 3) Mnożenie L_{i-2}(x) przez skalar -(i-1)^2
    # 4) Sumujemy współczynniki przy tych samych potęgach x z wyników operacji 1), 2), 3)
# Podstawienie L_k(x) do wzoru na P_N(x)
# P_N(x) = Σ_{k=0}^N c_k * (Σ_{j=0}^k l_{k,j} * x^j)
# Musimy zebrać wszystkie wyrazy, które mnożą się przez x^m, gdzie m = 0, 1, ..., N.
# Stąd wzór na współczynnik a_m to:
# a_m = Σ_{k=0}^N c_k * l_{k,m} dla m = 0, 1, ..., N

def poly_multiply_by_x(p: list[float]) -> list[float]:
    #Np: [a, b, c] (ax^2+bx+c) -> [a, b, c, 0] (ax^3+bx^2+cx+0)
    return p + [0.0]

def poly_multiply_by_scalar(p: list[float], scalar: float) -> list[float]:
    # Mnożenie wielomianu przez skalar
    return [coeff * scalar for coeff in p]

def poly_add(p1: list[float], p2: list[float]) -> list[float]:
    # Dodawanie dwóch wielomianów
    max_len = max(len(p1), len(p2))
    p1_extended = p1 + [0.0] * (max_len - len(p1))
    p2_extended = p2 + [0.0] * (max_len - len(p2))
    return [p1_extended[i] + p2_extended[i] for i in range(max_len)]

def calculate_polynomial_coefficients_from_laguerre(
        max_n: int) -> list[float]:

    # max_n - maksymalny stopień wielomianu, dla którego chcemy obliczyć współczynniki
    all_lk = []
    coeffs_l0 = [1.0]
    all_lk.append(coeffs_l0)
    if max_n == 0:
        return all_lk

    coeffs_l1 = [-1.0, 1.0]
    all_lk.append(coeffs_l1)
    if max_n == 1:
        return all_lk

    for i in range(2, max_n + 1):

        coeffs_l_prev = all_lk[i - 1]
        coeffs_l_prev_prev = all_lk[i - 2]
        k = i - 1

        # Operacja 1) Mnożenie L_{i-1}(x) przez x
        coeffs_l_prev_multiplied_by_x = poly_multiply_by_x(coeffs_l_prev)

        # Operacja 2) Mnożenie L_{i-1}(x) przez skalar, czyli -(2*k + 1)
        scalar_1 = -(2 * k + 1)
        coeffs_l_prev_multiplied_by_scalar = poly_multiply_by_scalar(coeffs_l_prev, scalar_1)

        # Operacja 3) Mnożenie L_{i-2}(x) przez skalar, czyli -k^2
        scalar_2 = -(k ** 2)
        coeffs_l_prev_prev_multiplied_by_scalar = poly_multiply_by_scalar(coeffs_l_prev_prev, scalar_2)

        # Operacja 4) Sumowanie współczynników przy tych samych potęgach x
        sum_1_2 = poly_add(
            coeffs_l_prev_multiplied_by_x,
            coeffs_l_prev_multiplied_by_scalar
        )

        sum_1_2_3 = poly_add(
            sum_1_2,
            coeffs_l_prev_prev_multiplied_by_scalar
        )

        # Dodajemy obliczony wielomian L_k(x) do listy
        all_lk.append(sum_1_2_3)

    return all_lk



# TODO: Użyć hornera do obliczenia wartości wielomianu w punkcie x.

