import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import inspect
from funzioni import *


def sostituisci_nan_con_vicino(vettore):
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
    odmr['freq'] = data[:, 0]  # in kHz
    odmr['ref']  = data[:, 1]   # in V
    odmr['od']   = data[:, 2]    # in V
    odmr['lock'] = data[:, 3]
    
    # Controlla e sostituisce i NaN in ciascun vettore
    odmr['freq'] = sostituisci_nan_con_vicino(odmr['freq'])
    odmr['ref'] = sostituisci_nan_con_vicino(odmr['ref'])
    odmr['od'] = sostituisci_nan_con_vicino(odmr['od'])
    odmr['lock'] = sostituisci_nan_con_vicino(odmr['lock'])

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
        p0[1] = 0.42  # sempre k | un 2π ogni 15 kHz -> 2π/15
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

color = ['blue', 'green', 'orange', 'purple']

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
    plt.xlabel('Frequenza [kHz]')
    plt.ylabel(yaxis)
    if title:
        plt.title(title+fr" | D={dist}mm")
    else:
        plt.title(fr"D=${dist}$mm | {yaxis} vs Frequenza")
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