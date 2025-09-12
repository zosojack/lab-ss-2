import numpy as np

# Funzioni di BASE  + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

def sin_func(x, A, k, phi):
        return A * np.sin(k * x + phi)
    
def lin_func(x, m, q):
        return m*x + q
    
def quad_func(x, a, b, c):
        return a*x*x + b*x + c
    
def cub_func(x, a, b, c, d):
    return a*x*x*x + b*x*x + c*x + d

def quar_func(x, a, b, c, d, e):
    return a*x*x*x*x + b*x*x*x + c*x*x + d*x + e

def quin_func(x, a, b, c, d, e, f):
    return a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x + f

def exa_func(x, a, b, c, d, e, f, g):
    return a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x*x + f*x + g

def ept_func(x, a, b, c, d, e, f, g, h):
    return a*x**7 + b*x**6 + c*x**5 + d*x**4 + e*x**3 + f*x**2 + g*x + h

def esp_func(x, tau):
    return np.exp(-x/tau)

# Funzioni per il fitting di curve ODMR + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
# retta + seno costante
def lin_plus_sin_func(x, A, k, phi, m, q):
    return lin_func(x, m, q) + sin_func(x, A, k, phi)

# retta + seno variabile
def lin_plus_var_sin_func(x, A, k, phi, m, q, p):
    return lin_func(x, m, q) + sin_func(np.pow(x,p), A, k, phi)
# stessa ma con logaritmo
def lin_plus_var_sin_func_log(x, A, k, phi, m, q):
    return lin_func(x, m, q) + sin_func(np.log(x), A, k, phi)

# retta con pendenza modulata sinusoide
def lin_mod_sin_func(x, A, k, phi, q):
    return x * sin_func(x, A, k, phi) + q

# quadratica + sinusoide
def quad_plus_sin_func(x, A, k, phi, a, b, c):
    return quad_func(x, a, b, c) + sin_func(x, A, k, phi)

# quadratica + sinusoide variabile
def quad_plus_var_sin_func(x, A, k, phi, a, b, c, p):
    return quad_func(x, a, b, c) + sin_func(np.pow(x,p), A, k, phi)
# stessa ma con logaritmo
def quad_plus_var_sin_func_log(x, A, k, phi, a, b, c):
    return quad_func(x, a, b, c) + sin_func(np.log(x), A, k, phi)

# cubica + sinusoide
def cub_plus_sin_func(x, A, k, phi, a, b, c, d):
    return cub_func(x, a, b, c, d) + sin_func(x, A, k, phi)

# cubica + sinusoide variabile
def cub_plus_var_sin_func(x, A, k, phi, a, b, c, d, p):
    return cub_func(x, a, b, c, d) + sin_func(x**p, A, k, phi)
# stessa ma con logaritmo
def cub_plus_var_sin_func_log(x, A, k, phi, a, b, c, d):
    return cub_func(x, a, b, c, d) + sin_func(np.log(x), A, k, phi)

# quartica + sinusoide
def quar_plus_sin_func(x, A, k, phi, a, b, c, d, e):
    return quar_func(x, a, b, c, d, e) + sin_func(x, A, k, phi)

# quintica + sinusoide
def quin_plus_sin_func(x, A, k, phi, a, b, c, d, e, f):
    return quin_func(x, a, b, c, d, e, f) + sin_func(x, A, k, phi)

# exa
def exa_plus_sin_func(x, A, k, phi, a, b, c, d, e, f, g):
    return exa_func(x, a, b, c, d, e, f, g) + sin_func(x, A, k, phi)

def ept_plus_sin_func(x, A, k, phi, a, b, c, d, e, f, g, h):
    return ept_func(x, a, b, c, d, e, f, g, h) + sin_func(x, A, k, phi)


# zio pera iniziano a diventare terribili

def lin_plus_var_sin_func_log_smorz(x, A, k, phi, a, b, tau):
    return esp_func(x, tau)* ( lin_func(x, a, b) + sin_func(np.log(x), A, k, phi) )

def quad_plus_var_sin_func_log_smorz(x, A, k, phi, a, b, c, tau):
    return esp_func(x, tau)* ( quad_func(x, a, b, c) + sin_func(np.log(x), A, k, phi) )

def cub_plus_var_sin_func_log_smorz(x, A, k, phi, a, b, c, d, tau):
    return esp_func(x, tau)* ( cub_func(x, a, b, c, d) + sin_func(np.log(x), A, k, phi) )


# prova nuova
def lin_plus_cos(x, A, k, phi, m, q):
    return lin_func(x, m, q) + A * np.cos(k / (x + phi))

def quad_plus_cos(x, A, k, phi, a, b, c):
    return quad_func(x, a, b, c) + A * np.cos(k / (x + phi))

def cub_plus_cos(x, A, k, phi, a, b, c, d):
    return cub_func(x, a, b, c, d) + A * np.cos(k / (x + phi))


# seno variabile
def sin_var_func(t, A, T0, k, dx, dy, alpha):
    """
    Seno variabile con smorzamento esponenziale.
    
    Args:
        t (float or np.ndarray): Tempo o array di tempi.
        A (float): Ampiezza del segnale.
        T0 (float): Periodo iniziale del segnale.
        k (float): Coefficiente di modulazione della frequenza.
        dy (float): Offset costante.
        dx (float): Offset costante.
        α (float): Coefficiente di smorzamento esponenziale.

    Returns:
        float or np.ndarray: Valore del segnale al tempo t.
    """
    return A * np.sin(2 * np.pi * (t-dx) / (T0 + k * (t-dx))) * np.exp(-alpha * (t-dx)) + dy

# uguale ma con seno quadro
def sin_quadro_var_func(t, A, T0, k, dx, dy, alpha):
    """
    Seno variabile con smorzamento esponenziale.
    
    Args:
        t (float or np.ndarray): Tempo o array di tempi.
        A (float): Ampiezza del segnale.
        T0 (float): Periodo iniziale del segnale.
        k (float): Coefficiente di modulazione della frequenza.
        dy (float): Offset costante.
        dx (float): Offset costante.
        α (float): Coefficiente di smorzamento esponenziale.

    Returns:
        float or np.ndarray: Valore del segnale al tempo t.
    """
    return A * (np.sin(2 * np.pi * (t-dx) / (T0 + k * (t-dx)))) * (np.sin(2 * np.pi * (t-dx) / (T0 + k * (t-dx)))) * np.exp(-alpha * (t-dx)) + dy