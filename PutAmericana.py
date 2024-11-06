import numpy as np
import pandas as pnd
from scipy.stats import norm

def put_americana_trinomial(S0, K, T, r, sigma, N):
    """
    Valora una opción put americana utilizando un árbol trinomial.
    
    Parámetros:
    S0 : float - Precio inicial del activo subyacente
    K : float - Precio de ejercicio de la opción
    T : float - Tiempo hasta el vencimiento en años
    r : float - Tasa de interés libre de riesgo
    sigma : float - Volatilidad del activo subyacente
    N : int - Número de pasos en el árbol
    
    Retorna:
    float - Valor de la opción put americana
    """
    # Parámetros del árbol
    dt = T / N  # Intervalo de tiempo por paso
    u = np.exp(sigma * np.sqrt(3 * dt))  # Factor de subida
    d = 1 / u  # Factor de bajada
    m = 1  # Factor neutro (sin movimiento)
    
    # Probabilidades de transición
    pu = 1/6
    pd = 1/6
    pm = 2/3
    
    # Inicialización del árbol de precios del activo subyacente
    precios_activo = np.zeros((2 * N + 1, N + 1))
    precios_activo[N, 0] = S0  # El precio inicial está en el centro del árbol
    
    for j in range(1, N + 1):
        for i in range(N - j, N + j + 1):
            precios_activo[i, j] = S0 * (u ** max(j - (N - i), 0)) * (d ** max((N + i) - j, 0))
    
    # Inicialización del árbol de valores de la opción
    valores_opcion = np.zeros((2 * N + 1, N + 1))
    
    # Valor de la opción en el último paso (valor intrínseco de la put)
    for i in range(2 * N + 1):
        valores_opcion[i, N] = max(K - precios_activo[i, N], 0)
    
    # Retropropagación para calcular el valor de la opción en t=0
    factor_descuento = np.exp(-r * dt)
    for j in range(N - 1, -1, -1):
        for i in range(N - j, N + j + 1):
            valor_mantener = (pu * valores_opcion[i - 1, j + 1] +
                              pm * valores_opcion[i, j + 1] +
                              pd * valores_opcion[i + 1, j + 1]) * factor_descuento
            valor_ejercicio = max(K - precios_activo[i, j], 0)
            valores_opcion[i, j] = max(valor_mantener, valor_ejercicio)
    
    return valores_opcion[N, 0]


# Ejemplo de uso

s0 = 30  # Precio inicial del activo
K = 35   # Strike
r = 0.05  # Tasa libre de riesgo (12%)
sigma = 0.2  # Volatilidad (20%)
T = 1  # Tiempo hasta la madurez (1 año)
N = 4  # Número de pasos en el árbol

# Valor de la opción put americana
poblacion = 20 

# Calculo el precio de accion para 20 valores de s0
lista_precios_tm = []
for i in range(poblacion):
    precio_put = put_americana_trinomial(s0+i,K, T, r, sigma, N)
    lista_precios_tm.append(precio_put)
print(f"Precios de opción put con trinomial normal: {lista_precios_tm}")


print("El valor de la opción put americana es:", lista_precios_tm)

# Blackscholes 
def black_scholes_put(S, K, T, r, sigma):

    # Cálculo de d1 y d2
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    
    # Cálculo del valor de la opción put
    put_price = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    
    return put_price

# Ejemplo de uso
s0 = 40  # Precio inicial del activo
K = 35   # Strike
r = 0.05  # Tasa libre de riesgo (12%)
sigma = 0.2  # Volatilidad (20%)
T = 1  # Tiempo hasta la madurez (1 año)

# Calculo el precio de una put con black scholes para 20 valores de s0
lista_precios_bs = []
for i in range(poblacion):
    put_price = black_scholes_put(s0+i, K, T, r, sigma)
    lista_precios_bs.append(put_price)

print("El valor de la opción put es:", lista_precios_bs)

# Calculo el error absoluto promedio de los 20 experimentos
error_absoluto_prom = 0
for i in range(len(lista_precios_tm)):
    error_absoluto_prom += abs(lista_precios_tm[i] - lista_precios_bs[i])

error_absoluto_prom = error_absoluto_prom/len(lista_precios_tm)
print(" Error absoluto promedio = ", error_absoluto_prom)
