import numpy as np
import pandas as pnd
from scipy.stats import norm

# Esta funcion calcula el precio de opcion 
# para los elementos de la malla (Es un arbol trinomial normal)
def calcular_elemento_malla(S0, k_pasos, dt, sigma):
    
    # El dt que se pasa aqui seria el valor de tiempo que pasa entre N-1 y N. Entonces ese tiempo se divide por los k_pasos
    dt_malla = dt / k_pasos
    # Parametros del arbol trinomial
    h = (sigma * np.sqrt(3*dt_malla))/2

    pu = 1/6
    pd = 1/6
    pm = 2/3

    # Inicializacion de los precios del activo en el arbol

    # S es un arreglo mas bonito
    S = np.zeros((2 * k_pasos + 1, k_pasos + 1 )) # Parametro es (Cantidad de filas, Cantidad de columnas)

    S[k_pasos, 0] = S0  # Nodo inicial para calcular

    # Rellenamos los precios del activo en cada nodo

    for t in range(1, k_pasos + 1):
        for i in range(k_pasos - t, k_pasos + t + 1, 1):
            if i >= k_pasos:
                S[i, t] = S[k_pasos, 0] - h * (i - k_pasos)
            elif i < k_pasos:
                S[i, t] = S[k_pasos, 0] + h * (k_pasos - i)
    
    # Inicializacion de los valores de la opcion en el tiempo de madurez
    
    precios_opcion = np.zeros((2 * k_pasos + 1, k_pasos + 1)) # Mismo tamaño que el precio, pero se recorre de adelante hacia atras

    # Precio opcion para una put (Con log aplicado, pues S0 es con log)
    precios_opcion[:, k_pasos] = np.maximum(K - np.exp(S[:, k_pasos]), 0)
    #print(precios_opcion)
    # Retropropagacion de los valores de la opcion en el arbol trinomial
    for t in range(k_pasos - 1, -1, -1):
        for i in range(k_pasos - t, k_pasos + t + 1, 1):
            #print(k_pasos)
            valor_de_continuacion = pu * precios_opcion[i+1, t+1] + pm * precios_opcion[i, t+1] + pd * precios_opcion[i-1, t+1]
            valor_de_continuacion *= np.exp(-r * T)

            precios_opcion[i, t] = np.maximum(K - np.exp(S[i, t]), valor_de_continuacion)
    
    #print(precios_opcion)
    return precios_opcion[k_pasos, 0]


def trinomial_grande_put(S0, N, dt, sigma):
    # Parametro utilizado en el paper
    h = sigma * np.sqrt(3*dt)
    
    # Probabilidades en el arbol trinomial

    pu = 1/6
    pd = 1/6
    pm = 2/3

    # Inicializacion de los precios del activo en el arbol
    S = np.zeros((2*N + 1, N + 1)) # Parametro es (Cantidad de filas, Cantidad de columnas)
    # Valor del nodo S0, aplicandole el logaritmo
    S[N, 0] = np.log(S0)

    # Rellenamos los precios del activo en cada nodo. 
    for t in range(1, N+1):
        for i in range(N - t, N + t + 1, 1):
            if i >= N:
                # Ver que N - i da positivo, y se suma la cantaidad de h's ( Sube la variable X*)
                S[i, t] = S[N, 0] + h * (i-N)
            elif i < N:
                #  Ver que N - i da positivo, y se resta la cantidad de h's ( Baja la variable X*)
                S[i, t] = S[N, 0] - h * (N-i)

    # Pasos de la malla
    k = 4 

    # Calculamos los nodos en los que hay que hacer la malla.

    nodos_malla_fina = []
    nodo_superior_k = 0
    nodo_inferior_k = 0
    # Aqui se eligen los nodos en los que se realizara la malla
    for i in range(1, N*2):
        
        # Chequeo entre que nodos de precio se encuentra el strike K
        if np.exp(S[i, N-1]) <= K and np.exp(S[i + 1, N-1]) >= K:
            nodo_superior_k = S[i+1, N-1]
            nodo_inferior_k = S[i, N-1]
            
            # Agrego los nodos que encierran a K de manera ordenada.
            if i > 0:
                nodos_malla_fina.append(i-1)
            nodos_malla_fina.append(i)
            nodos_malla_fina.append(i+1)
            if i+2 <= N+1:
                nodos_malla_fina.append(i+2)
    
    # Aca obtengo los precios de la malla fina, para usarlos como S0 en calcular_elemento_malla
    precios_nodos_malla_fina = [S[i, N-1] for i in nodos_malla_fina]

    # Calculo los valores de opcion del arbol grande, de los nodos de la malla fina
    mallas = []
    for i in range(len(precios_nodos_malla_fina)):
        #
        valor_inicial = precios_nodos_malla_fina[i]
        mallas.append(calcular_elemento_malla(valor_inicial, k, dt, sigma))
    
    # Inicializacion de los valores de la opcion en el tiempo de la madurez

    opcion_precios = np.zeros((2*N+1, N+1))
    opcion_precios[:, N] = np.maximum(K - np.exp(S[:, N]), 0)

    # Retropropagacion de los valores de la opcion en el arbol trinomial

    for t in range(N-1, -1, -1):
        for i in range(N-t, N+t+1, 1):
            valor_de_continuacion = pu * opcion_precios[i+1, t+1] + pm * opcion_precios[i, t+1] + pd * opcion_precios[i-1, t+1]
            valor_de_continuacion *= np.exp(-r * T)

            if (i in nodos_malla_fina) and t == N-1:
                # Encuentro el indice del precio de opcion y lo asigno
                aux = nodos_malla_fina.index(i)
                opcion_precios[i, t] = mallas[aux] 
            else:
                opcion_precios[i, t] = np.maximum(K - np.exp(S[i, t]), valor_de_continuacion)
    
    # Imprimir el array de opción usando pandas para un mejor formato
    # option_df = pnd.DataFrame(opcion_precios)
    # print("Precio opciones")
    # print(option_df)
    # price_df = pnd.DataFrame(S)
    # print("Arbol Precios ")
    # print(price_df)
    return opcion_precios[N, 0]


# Ejemplo de uso
s0 = 30  # Precio inicial del activo
K = 35   # Strike
r = 0.05  # Tasa libre de riesgo (12%)
sigma = 0.2  # Volatilidad (20%)
T = 1  # Tiempo hasta la madurez (1 año)
N = 4  # Número de pasos en el árbol

# Cantidad de experimentos
poblacion = 20

# Paso de tiempo
dt = T / N

# Precio de una opción put americana utilizando el árbol trinomial adaptativo
lista_precios_tm = []
for i in range(poblacion):
    precio_put = trinomial_grande_put(s0+i,N,dt,sigma)
    lista_precios_tm.append(precio_put)

print(f"Precios de opción put con trinomial con malla: {lista_precios_tm}")

# Put black scholes
def black_scholes_put(S, K, T, r, sigma):

    # Cálculo de d1 y d2
    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)
    
    # Cálculo del valor de la opción put
    put_price = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    
    return put_price

# Precios de opcion put americana usando black scholes
lista_precios_bs = []
for i in range(poblacion):
    put_price = black_scholes_put(s0+i, K, T, r, sigma)
    lista_precios_bs.append(put_price)
print("Precios con black scholes", lista_precios_bs)

# Calculo del valor absoluto
error_absoluto_prom = 0
for i in range(len(lista_precios_tm)):
    error_absoluto_prom += abs(lista_precios_tm[i] - lista_precios_bs[i])

error_absoluto_prom = error_absoluto_prom/len(lista_precios_tm)
print(" Error absoluto promedio = ", error_absoluto_prom)
    