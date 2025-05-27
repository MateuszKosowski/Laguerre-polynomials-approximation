import numpy as np

# Wzór wielomianu aproksymacyjnego Laguerre'a:
# L_n(x) = sum(c_k * L_k(x)), gdzie c_k to współczynniki aproksymacji

def laguerre_polynomial(n: int, x: float) -> float:
    # n - stopień wielomianu
    # x - punkt, w którym obliczamy wielomian
    # Zwracamy wartość wielomianu Laguerre'a dla danego n i x


    if n == 0:
        return 1.0
    elif n == 1:
        return 1.0 - x
    else:
        leguerre_prev_prev = 1.0  # Wielomian stopnia 0
        leguerre_prev = 1.0 - x  #  Wielomian stopnia 1

        for i in range(2, n + 1):

            # Iteracyjne obliczanie wielomianu Laguerre'a n stopnia za pomocą wzoru rekurencyjnego
            # Wzór rekurencyjny:
            # L_n(x) = ((2n-1-x)L_{n-1}(x) - (n-1)L_{n-2}(x)) / n
            leguerre_current = ((2 * i - 1 - x) * leguerre_prev - (i - 1) * leguerre_prev_prev) / i

            # Aktualizujemy poprzednie wielomiany
            leguerre_prev_prev, leguerre_prev = leguerre_prev, leguerre_current

        return leguerre_prev



