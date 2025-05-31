import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Zaimportuj funkcje z Twojego pliku leguerre.py
# Upewnij się, że leguerre.py jest w tym samym katalogu
from leguerre import (
    calculate_laguerre_polynomial_value,  # Choć bezpośrednio nie użyjemy jej w UI, jest potrzebna dla c_coeffs i reszty
    calculate_c_coeffs,
    calculate_polynomial_coefficients_from_laguerre,
    create_coefficients_of_power_polynomial,
    calculate_approximating_polynomial_value_by_horner
)


# --- Definicje Funkcji do Aproksymacji ---
def linear_func(x, a=1, b=0):
    return a * x + b


def abs_func(x):
    # Dla wielomianów Laguerre'a, które są zdefiniowane na [0, +nieskończoność),
    # |x| na tym przedziale to po prostu x.
    if x < 0:
        return -x  # Ogólna definicja, ale Laguerre jest dla x >= 0
    return x


def poly_func_example(x, a=1, b=-4, c=2):  # Przykład: L_2^*(x) = x^2 - 4x + 2
    return a * x ** 2 + b * x + c


def trig_func_sin(x):
    return np.sin(x)


def trig_func_cos(x):
    return np.cos(x)


def exp_func_simple(x):
    return np.exp(-x / 2)  # Funkcja dobrze zachowująca się dla Laguerre


def custom_L0(x):  # L_0^*(x)
    return calculate_laguerre_polynomial_value(0, x)


def custom_L1(x):  # L_1^*(x)
    return calculate_laguerre_polynomial_value(1, x)


def custom_L2(x):  # L_2^*(x)
    return calculate_laguerre_polynomial_value(2, x)


# Słownik funkcji dostępnych dla użytkownika
AVAILABLE_FUNCTIONS = {
    "Liniowa: 2x + 1": lambda x: linear_func(x, a=2, b=1),
    "Liniowa: x": lambda x: linear_func(x, a=1, b=0),  # Dla |x| na [0, inf)
    "Wielomian: x^2 - 4x + 2 (L_2*)": lambda x: poly_func_example(x, a=1, b=-4, c=2),
    "Wielomian: 0.5x^3 - 2x^2 + x - 1": lambda x: 0.5 * x ** 3 - 2 * x ** 2 + x - 1,
    "Trygonometryczna: sin(x)": trig_func_sin,
    "Trygonometryczna: cos(x)": trig_func_cos,
    "Trygonometryczna złożona: sin(x) + 2": lambda x: np.sin(x) + 2,
    "Wykładnicza: exp(-x/2)": exp_func_simple,
    "Test: L_0*(x) = 1": custom_L0,
    "Test: L_1*(x) = x - 1": custom_L1,
    "Test: L_2*(x) = x^2 - 4x + 2": custom_L2,  # To samo co poly_func_example
}


# --- Funkcja Pomocnicza do Tworzenia Stringa Wielomianu (LaTeX) ---
def get_polynomial_latex_string(coeffs_a: list[float], precision: int = 4) -> str:
    if not coeffs_a:
        return "0"

    poly_str_parts = []
    # Iterujemy od a_N do a_0 (najwyższa potęga do najniższej)
    # coeffs_a to [a0, a1, ..., aN]
    N = len(coeffs_a) - 1

    for i in range(N, -1, -1):  # i to aktualna potęga x
        a_val = coeffs_a[i]

        # Sprawdź, czy sformatowany współczynnik jest numerycznie zerem
        # np. jeśli a_val = 0.00003, a precision = 4, to formatted_val_for_check = "0.0000"
        # i float(formatted_val_for_check) będzie 0.0
        formatted_val_for_check = f"{a_val:.{precision}f}"  # Używamy oryginalnej wartości do formatowania

        if float(formatted_val_for_check) == 0.0:
            # Jeśli sformatowany współczynnik jest zerem:
            # - Jeśli jest to wyraz wolny (i=0) i nic wcześniej nie zostało dodane do wielomianu
            #   (czyli wszystkie wyższe potęgi też były zerowe po sformatowaniu),
            #   to dodaj "0.0000" (lub odpowiednio sformatowane zero).
            # - W przeciwnym razie (wyższe potęgi lub są już inne termy), pomiń ten zerowy wyraz.
            if i == 0 and not poly_str_parts:
                poly_str_parts.append(f"{0.0:.{precision}f}")  # Dodaj sformatowane zero jako jedyny term
            continue  # Pomiń ten wyraz

        term = ""
        # Znak
        if a_val > 0 and poly_str_parts:  # Jeśli dodatni i nie jest to pierwszy term dodawany do listy
            term += "+ "
        elif a_val < 0:  # Jeśli ujemny (niezależnie czy pierwszy, czy kolejny)
            term += "- "
        # Jeśli a_val > 0 i jest to pierwszy term (poly_str_parts jest puste), znaku "+" nie dodajemy na początku

        # Wartość współczynnika (bezwzględna)
        abs_val_formatted_str = f"{abs(a_val):.{precision}f}"

        # Czy pokazać wartość współczynnika "1.0000"?
        # Nie, jeśli potęga i > 0 (czyli dla x, x^2, ...) ORAZ sformatowana wartość to "1.0000"
        # Tak, jeśli potęga i = 0 (wyraz wolny)
        if i > 0 and abs_val_formatted_str == f"{1.0:.{precision}f}":
            # Nie dodawaj "1.0000" dla x, x^2, itp.
            # Jeśli `term` jest pusty (pierwszy wyraz, dodatni, np. x^k), to dobrze.
            # Jeśli `term` ma już znak ("- "), to też dobrze (np. -x^k).
            pass
        else:
            term += abs_val_formatted_str

        # Zmienna x i potęga
        if i > 0:  # Dla x, x^2, ...
            # Jeśli term jest pusty (np. pierwszy wyraz to x^k)
            # LUB jeśli term kończy się znakiem (np. "- " dla -x^k)
            # to nie dodawaj spacji przed 'x'.
            # LaTeX sam sobie poradzi: "1.23x" jest ok, "x" jest ok, "-x" jest ok, "- 1.23x" jest ok.
            term += "x"
            if i > 1:  # Dla x^2, x^3, ...
                term += f"^{{{i}}}"

        poly_str_parts.append(term)

    if not poly_str_parts:
        # Jeśli po całej pętli nic nie zostało dodane (wszystkie współczynniki były zerowe po sformatowaniu)
        return "0"  # Zwróć prosty string "0"

    final_polynomial_string = " ".join(poly_str_parts)

    # Usuń wiodący plus, jeśli jest (np. gdy pierwszy term był dodatni i nie był to `a_N`)
    # To się nie powinno zdarzyć przy obecnej logice budowania `term`,
    # bo `+ ` jest dodawane tylko jeśli `poly_str_parts` nie jest puste.
    # Ale dla pewności:
    if final_polynomial_string.startswith("+ "):
        final_polynomial_string = final_polynomial_string[2:]

    return final_polynomial_string if final_polynomial_string else "0"  # Ostateczne zabezpieczenie


# --- Główna część aplikacji Streamlit ---
st.set_page_config(layout="wide")
st.title("Aproksymacja Funkcji Wielomianami Laguerre'a")

st.sidebar.header("Ustawienia Aproksymacji")

# Wybór funkcji
selected_func_name = st.sidebar.selectbox(
    "Wybierz funkcję do aproksymacji:",
    list(AVAILABLE_FUNCTIONS.keys())
)
func_to_approximate = AVAILABLE_FUNCTIONS[selected_func_name]

# Parametry aproksymacji
n_degree_approximation = st.sidebar.slider(
    "Stopień wielomianu aproksymującego (N):",
    min_value=0, max_value=20, value=5, step=1
)
num_gauss_nodes_for_c_k = st.sidebar.slider(
    "Liczba węzłów kwadratury Gaussa-Laguerre'a (dla c_k):",
    min_value=2, max_value=100, value=30, step=1
)

# Parametry wykresu
st.sidebar.header("Ustawienia Wykresu")
plot_x_max = st.sidebar.slider(
    "Górna granica osi X dla wykresu:",
    min_value=1.0, max_value=50.0, value=10.0, step=0.5
)
plot_points = st.sidebar.slider(
    "Liczba punktów na wykresie:",
    min_value=50, max_value=1000, value=400, step=50
)

# Wyświetlanie informacji o wyborze
st.header("Parametry")
st.markdown(f"""
- **Aproksymowana funkcja f(x):** `{selected_func_name}`
- **Stopień wielomianu aproksymującego N:** `{n_degree_approximation}`
- **Liczba węzłów Gaussa-Laguerre'a dla c_k:** `{num_gauss_nodes_for_c_k}`
""")

# Obliczenia
if st.button("Rozpocznij aproksymację"):
    st.header("Wyniki Aproksymacji")

    with st.spinner("Obliczanie współczynników c_k..."):
        try:
            calculated_c_coefficients = calculate_c_coeffs(
                func=func_to_approximate,
                n=n_degree_approximation,
                num_nodes_gauss=num_gauss_nodes_for_c_k
            )
        except Exception as e:
            st.error(f"Błąd podczas obliczania współczynników c_k: {e}")
            st.stop()

    with st.spinner("Generowanie współczynników L_k(x) w bazie potęgowej..."):
        try:
            all_lk_coeffs = calculate_polynomial_coefficients_from_laguerre(n_degree_approximation)
        except Exception as e:
            st.error(f"Błąd podczas generowania współczynników L_k(x): {e}")
            st.stop()

    with st.spinner("Transformacja współczynników c_k na a_m (baza potęgowa)..."):
        try:
            a_m_power_coefficients = create_coefficients_of_power_polynomial(
                calculated_c_coefficients,
                all_lk_coeffs
            )
        except Exception as e:
            st.error(f"Błąd podczas transformacji c_k na a_m: {e}")
            st.stop()

    # Wyświetlanie współczynników
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Obliczone współczynniki c_k")
        c_k_data = {"k": list(range(len(calculated_c_coefficients))),
                    "c_k": [f"{c:.6e}" for c in calculated_c_coefficients]}
        st.dataframe(c_k_data, height=min(300, 35 * (n_degree_approximation + 2)))

    with col2:
        st.subheader("Współczynniki a_m (baza potęgowa)")
        a_m_data = {"m": list(range(len(a_m_power_coefficients))), "a_m": [f"{a:.6e}" for a in a_m_power_coefficients]}
        st.dataframe(a_m_data, height=min(300, 35 * (n_degree_approximation + 2)))

    st.subheader(f"Wielomian aproksymujący P_{n_degree_approximation}(x)")
    poly_latex_str = get_polynomial_latex_string(a_m_power_coefficients)
    st.latex(poly_latex_str)

    # Przygotowanie danych do wykresu
    x_vals_plot = np.linspace(0, plot_x_max, plot_points)  # Wielomiany Laguerre'a są dla x >= 0

    y_original = np.zeros_like(x_vals_plot)
    for i, x_val in enumerate(x_vals_plot):
        try:
            y_original[i] = func_to_approximate(x_val)
        except Exception:  # np.log(0)
            y_original[i] = np.nan

    y_approximated_horner = np.array([
        calculate_approximating_polynomial_value_by_horner(x, a_m_power_coefficients)
        for x in x_vals_plot
    ])

    # Wykres
    st.subheader("Wykres funkcji oryginalnej i aproksymacji")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x_vals_plot, y_original, label=f"f(x) = {selected_func_name}", color="blue")
    ax.plot(x_vals_plot, y_approximated_horner, label=f"P_{n_degree_approximation}(x) (Horner)", linestyle="--",
            color="red")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Porównanie funkcji oryginalnej i aproksymacji Laguerre'a")
    ax.legend()
    ax.grid(True)
    ax.set_ylim(bottom=min(np.nanmin(y_original) if not np.all(np.isnan(y_original)) else 0,
                           np.nanmin(y_approximated_horner) if not np.all(np.isnan(y_approximated_horner)) else 0) - 1,
                top=max(np.nanmax(y_original) if not np.all(np.isnan(y_original)) else 0,
                        np.nanmax(y_approximated_horner) if not np.all(
                            np.isnan(y_approximated_horner)) else 0) + 1)  # Dynamiczne Y lim
    st.pyplot(fig)

    # Obliczanie błędu
    st.subheader("Ocena Błędu Aproksymacji")
    # Usuń NaN wartości, które mogły powstać w y_original (np. log(0))
    valid_indices = ~np.isnan(y_original) & ~np.isnan(y_approximated_horner)
    if np.any(valid_indices):
        y_orig_valid = y_original[valid_indices]
        y_approx_valid = y_approximated_horner[valid_indices]

        errors = np.abs(y_orig_valid - y_approx_valid)
        max_abs_error = np.max(errors)
        mean_squared_error = np.mean(errors ** 2)

        st.metric(label="Maksymalny Błąd Absolutny (MAE)", value=f"{max_abs_error:.4e}")
        st.metric(label="Błąd Średniokwadratowy (MSE)", value=f"{mean_squared_error:.4e}")


    else:
        st.warning("Nie można obliczyć błędu - brak poprawnych punktów danych.")

else:
    st.info("Kliknij przycisk 'Rozpocznij aproksymację', aby zobaczyć wyniki.")

st.sidebar.markdown("---")
st.sidebar.markdown("Wariant 3: Wielomiany Laguerre'a.")