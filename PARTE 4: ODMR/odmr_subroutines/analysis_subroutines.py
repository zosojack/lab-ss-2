import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import inspect
from odmr_subroutines.funzioni import *

# questa funzione è utile al momento della lettura dei file
def _sostituisci_nan_con_vicino(vettore):
    """
    Sostituisce i valori NaN in un array con il valore a sinistra.
    Se il primo elemento è NaN, usa il valore a destra.
    """
    for i in range(len(vettore)):
        if np.isnan(vettore[i]):
            if i == 0:  # Se è il primo elemento, usa il valore a destra
                vettore[i] = vettore[i + 1]
            else:  # Altrimenti usa il valore a sinistra
                vettore[i] = vettore[i - 1]
    return vettore

def leggi_file_odmr(nome):
    odmr = {
        'freq': np.ndarray([]),
        'ref':  np.ndarray([]),
        'od':   np.ndarray([]),
        'lock': np.ndarray([])
    }
    
    data = np.loadtxt(nome, skiprows=2, delimiter='|')
    odmr['freq'] = data[:, 0]  # in MHz
    odmr['ref']  = data[:, 1]   # in V
    odmr['od']   = data[:, 2]    # in V
    odmr['lock'] = data[:, 3]
    
    # Controlla e sostituisce i NaN in ciascun vettore
    odmr['freq'] = _sostituisci_nan_con_vicino(odmr['freq'])
    odmr['ref'] = _sostituisci_nan_con_vicino(odmr['ref'])
    odmr['od'] = _sostituisci_nan_con_vicino(odmr['od'])
    odmr['lock'] = _sostituisci_nan_con_vicino(odmr['lock'])

    return odmr

# lo spettro lock-in non ha il fondo attorno a 0, occorre una traslazione verticale
def trasla_spettro_lock_in(odmr):
    """
    Trasla verticalmente il segnale lock-in della sua stessa media.
    """
    media = np.mean(odmr['lock'])
    odmr['lock'] -= media
    return odmr

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
    
def fit_sin_odmr(odmr, key: str, region: tuple = None, p0 = None, sinp0=None, fit_func = lin_plus_var_sin_func):
    
    if region:
        mask = (odmr['freq'] >= region[0]) & (odmr['freq'] <= region[1])
        freq = odmr['freq'][mask].copy()
        key_data = odmr[key][mask].copy()
    else:
        freq = odmr['freq']
        key_data = odmr[key].copy()
    
    
    print("FUNZIONE DI FIT:", fit_func.__name__)
    
    # Get the number of parameters for the fit function
    num_params = len(inspect.signature(fit_func).parameters)
    
    # bounds
    if fit_func == lin_plus_var_sin_func or fit_func == quad_plus_var_sin_func or fit_func == cub_plus_var_sin_func:
        upp = [1, np.inf, np.pi] + [np.inf]*(num_params-4-1) + [1]
        low = [-1, -np.inf, -np.pi] + [-np.inf]*(num_params-4-1) + [0]
    elif fit_func == lin_plus_var_sin_func_log_smorz or fit_func == quad_plus_var_sin_func_log_smorz or fit_func == cub_plus_var_sin_func_log_smorz:
        upp = [1, np.inf, np.pi] + [np.inf]*(num_params-4)
        low = [-1, -np.inf, -np.pi] + [-np.inf]*(num_params-4)
    else:
        upp = [np.inf] * (num_params-1)
        low = [-np.inf] * (num_params-1)
    
    
    if p0 is None:
        # Generate initial p0 with appropriate number of values
        p0 = [0.1] * (num_params-1)
        p0[0] = 0.005  # sempre a
        p0[1] = 0.42  # sempre k | un 2π ogni 15 MHz -> 2π/15
        p0[2] = 0  # sempre phi
        
        
    # initial parameters
    if sinp0 is not None:
        p0 = [0.1] * (num_params-1)
        p0[0] = sinp0[0]
        p0[1] = sinp0[1]
        p0[2] = sinp0[2]
        # fisso i parametri iniziali
        upp[0], upp[1], upp[2] = p0[0] + 1e-10, p0[1] + 1e-10, p0[2] + 1e-10
        low[0], low[1], low[2] = p0[0] - 1e-10, p0[1] - 1e-10, p0[2] - 1e-10
        
        
    if fit_func == lin_plus_var_sin_func or fit_func == quad_plus_var_sin_func or fit_func == cub_plus_var_sin_func:
        p0[-1] = 0.9
    elif fit_func == lin_plus_var_sin_func_log_smorz or fit_func == quad_plus_var_sin_func_log_smorz or fit_func == cub_plus_var_sin_func_log_smorz:
        p0[-1] = 1000
    
    popt, _ = curve_fit(fit_func, freq, key_data, p0=p0, maxfev=500000, bounds=(low,upp))
    
    fit = fit_func(freq, *popt)
    
    print("PARAMETRI FIT:", popt)
    
    return popt, (freq, fit)


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +


# PLOTTER #

color = ['blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan', 'magenta', 'yellow', 'teal', 'navy', 'maroon', 'lime', 'coral', 'gold', 'silver']

def plot_odmr(odmr, key: str, title=None, dist=0, fit_curves=None):
    
    if key == 'lock':
        yaxis = 'Segnale Lock-In [V]'
    else:
        yaxis = 'Intensità PL [V]'
        
    if dist == 1000:
        dist = '\infty'
    
    plt.figure(figsize=(10, 3), dpi=200)
    plt.plot(odmr['freq'], odmr[key], color='red')    
    if fit_curves is not None:
        for i, fit_curve in enumerate(fit_curves):
            plt.plot(fit_curve[0], fit_curve[1], color=color[i], linestyle='--', label=f'Fit {i+1}')
            plt.legend()
    plt.xlabel('Frequenza [MHz]')
    plt.ylabel(yaxis)
    if title:
        plt.title(title+fr" | D={dist}mm")
    else:
        plt.title(rf"D=${dist}$mm | {yaxis} vs Frequenza")
    plt.grid(linestyle='--')
    plt.tight_layout()
    plt.show()
    
    
def remove_background(odmr, key: str, fit_func, fit_curve, popt, regione=None):
    """
    Rimuove il background da un dataset ODMR.
    """
    if regione:
        mask = (odmr['freq'] >= regione[0]) & (odmr['freq'] <= regione[1])
    else:
        mask = np.arange(len(odmr['freq']))
        
    curva_da_sottrarre = fit_func(odmr['freq'][mask], *popt)
    odmr[key][mask] -= curva_da_sottrarre
    
    
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

# ODMR #
def lorentziana(nu, A, gamma, nu0, q):
    """
    Funzione di Lorentziana più background lineare.
    A: ampiezza della Lorentziana
    gamma: half-width at half-maximum (HWHM)
    nu0: centro della Lorentziana
    nu: frequenza (variabile indipendente)
    """
    return -A / (1 + ((nu - nu0)/gamma)**2) + q

def lorentziana_più_quad(nu, A, gamma, nu0, a, b, c):
    """
    Funzione di Lorentziana più background lineare.
    A: ampiezza della Lorentziana
    gamma: half-width at half-maximum (HWHM)
    nu0: centro della Lorentziana
    a: coefficiente angolare della retta di background
    b: coefficiente angolare della quadratica di background
    c: intercetta della quadratica di background
    nu: frequenza (variabile indipendente)
    """
    return -A / (1 + ((nu - nu0)/gamma)**2) + a * nu**2 + b * nu + c

def fitta_deep_separati(odmr, N=1, regione=None, nu0s=None):
    """
    Esegue il fit della risonanza profonda (deep) in uno spettro ODMR.
    """
    
    params = []
    curva = np.zeros_like(odmr['freq'])
    n_curve = []
    
    print("Fit separato:")
    for i in range(0, N):
        singola_lorentziana = []

        mask = (odmr['freq'] >= nu0s[i]-20) & (odmr['freq'] <= nu0s[i]+20)
        freq = odmr['freq'][mask].copy()
        od_data = odmr['od'][mask].copy()
        
        help = np.zeros_like(freq)
        
        # SE LORENTZ + QUADRATICA
        '''
        if nu0s is None:
            p0 = [0.5, 5, 2870, 0., 0., 0.]  # A, gamma, nu0, a, b, c
            low = [0, 0, 2860, -np.inf, -np.inf, -1]
            upp = [1, 20, 2880, np.inf, np.inf, 1]
        else:
            p0 = [0.5, 5, nu0s[i], 0., 0., 0.]  # A, gamma, nu0, a, b, c
            low = [0, 0, nu0s[i]-5, -np.inf, -np.inf, -1]
            upp = [1, 10, nu0s[i]+5, np.inf, np.inf, 1]
        '''
        # SENZA PARABOLA
        if nu0s is None:
            p0 = [0.5, 5, 2870, 0., 0.]  # A, gamma, nu0, m, q
            low = [0, 0, 2860, -np.inf, -1]
            upp = [1., 30, 2880, np.inf, 1]
        else:
            p0 = [0.5, 5, nu0s[i], 0.]  # A, gamma, nu0, q
            low = [0, 0, nu0s[i]-20, -1]
            upp = [1., 30, nu0s[i]+20, 1] 

        popt, pcov = curve_fit(lorentziana, freq, od_data, p0=p0, maxfev=500000, bounds=(low,upp))

        singola_lorentziana.append( (popt, pcov) )
        curva[mask] += lorentziana(freq, *popt)
        help = lorentziana(freq, *popt)
        n_curve.append((freq,help))
        
        print(f"{i}. A={popt[0]}, gamma={popt[1]}, nu0={popt[2]}")

        params.append(singola_lorentziana)
    
    # così ritorno una singola curva
    #return params, (odmr['freq'], curva)
    
    # così ritorno N curve
    return params, n_curve

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def fitta_deep_insieme(odmr, N=1, regione=None, nu0s=None):
    """
    Esegue il fit della risonanza profonda (deep) in uno spettro ODMR.
    """

    def n_lorentziane(nu, *params):
        """
        Somma di N lorentziane più traslazione.
        params: lista di parametri [A1, gamma1, nu01, A2, gamma2, nu02, ..., q]
        """
        result = np.zeros_like(nu)
        for i in range(N):
            A = params[3*i]
            gamma = params[3*i + 1]
            nu0 = params[3*i + 2]
            result += -A / (1 + ((nu - nu0)/gamma)**2)
        # Traslazione
        q = params[-1]
        result += q
        return result

    
    if regione:
        mask = (odmr['freq'] >= regione[0]) & (odmr['freq'] <= regione[1])
        freq = odmr['freq'][mask].copy()
        od_data = odmr['od'][mask].copy()
    else:
        freq = odmr['freq']
        od_data = odmr['od'].copy()
    
    help = np.zeros_like(freq)
    
    p0 = []
    low = []
    upp = []
    
    # parametri lorentziane
    for nu0 in nu0s:
        p0 += [0.5, 5, nu0]  # A, gamma, nu0
        low += [0, 0, nu0-25]
        upp += [10., 30, nu0+25]
    # traslazione
    p0 += [0.]
    low += [-1]
    upp += [+1]

    popt, pcov = curve_fit(n_lorentziane, freq, od_data, p0=p0, maxfev=500000, bounds=(low,upp))

    print("Fit insieme:")
    for i in range(0, N):
        print(f"{i}. A={popt[3*i]}, gamma={popt[3*i+1]}, nu0={popt[3*i+2]}")
    
    # così ritorno N curve
    return (popt, pcov), (freq, n_lorentziane(freq, *popt))