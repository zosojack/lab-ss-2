import numpy as np
import matplotlib.pyplot as plt

# individua la regione dello sweep
def trova_regione (ch_B):
    # CONSIDERAZIONI
    # il ch_B ha fa un saltello a 2V
    # dopodiché oscilla attorno a un valore costante (1V)
    # quindi risaltella a 2 quando finisce
    # noi vogliamo la regione contenuta tra questi due saltelli

    # PIANO D'AZIONE
    # prendo il massimo all'inizio e il massimo alla fine con argmax
    # seleziono la mask tra questi due
    # poi faccio una mask dei valori sotto a 1.5V dentro a questa regione
    
    prima_metà = np.argmax(ch_B[:len(ch_B)//2])  # massimo nella prima metà
    seconda_metà = np.argmax(ch_B[len(ch_B)//2:]) + len(ch_B)//2  # massimo nella seconda metà
    mask = np.arange(prima_metà, seconda_metà)
    final_mask = mask[ch_B[mask] < 1.5]  # seleziona i valori sotto a 1.5V
    
    return final_mask

# ritaglia le intensità PL e converti tempi in frequenze nella regione
def intensità_e_frequenze_sweep (time, ch_A, ch_B):
    
    # io so che le frequenze hanno un passo di 0.25Hz, 
    # partono da 2500 e vanno a 3200
    # quindi posso impostare le frequenze a partire dal tempo
    # mi aspetto che ci siano 700*4 punti
    regione = trova_regione(ch_B)
    
    frequenze = [2500 + 0.25 * i for i in range(len(time[regione]))]  # 2500 + 0.25*i per i in [0, 2800] (len regione è il numero di punti [2800], non il delta T)
    intensità = ch_A[regione]  # intensità PL nella regione selezionata
    
    return frequenze, intensità


def trova_corrispondenza (intensità_ref, intensità_od, mean):
    # CONSIDERAZIONI
    # cerca di allineare le intensità ref con quelle odmr
    # e trova l'indice della frequenza corrispondente all'odmr
    
    # PIANO D'AZIONE
    # calcola la somma quadratica delle differenze tra le intensità
    # e trova il match con la somma minima
    
    # poiché l_ref è ritagliato, l_od è più lungo
    l_ref = len(intensità_ref)
    l_od = len(intensità_od)

    inizio_regione = 0 # quello da salvare
    fine_regione = l_ref - 1
    valore_precedente = +1E10
    
    # questo dice quante somme quadratiche vengono calcolate
    # e mi trasla i valori di ref su quelli di od
    for i in range(l_od - l_ref):
        somma_quadratica = 0
        # J gira effettivamente sui valori per calcolare le differenze

        if mean:
            # METODO ALTERNATIVO: 'REBIN' DIFFERENZA N VALORI ALLA VOLTA 
            for j in range(0, l_ref-5, 5):
                media_n_val_ref = np.mean(intensità_ref[j:j+5])
                media_n_val_od  = np.mean(intensità_od[j+i:j+5+i])
                somma_quadratica += (media_n_val_ref-media_n_val_od)**2
        else:
            # METODO BASICO: DIFFERENZA VALORE PER VALORE
            for j in range(l_ref):
                somma_quadratica += (intensità_ref[j]-intensità_od[j+i])**2
                
        if somma_quadratica < valore_precedente:
            inizio_regione = i
            fine_regione = i + l_ref
            
    return inizio_regione, fine_regione
    
    
def leggi_files (reference: str, odmr: str, mean):
    
    # leggo il reference
    data = np.loadtxt(reference, skiprows=9, delimiter=',')
    time_ref = data[:, 0]  # in s
    ch_A_ref = data[:, 1]  # in V - PL intensity
    ch_B_ref = data[:, 2]  # in V - reference

    # trovo la regione dello sweep e la ritaglio
    frequenze, intensità_ref = intensità_e_frequenze_sweep(time_ref, ch_A_ref, ch_B_ref)
    
    # leggo odmr
    data2 = np.loadtxt(odmr, skiprows=9, delimiter=',')
    ch_A_od = data2[:, 1]  # in V - PL intensity
    ch_B_od = data2[:, 2]  # in V - lock-in
    
    # cerco corrispondenzza tra intensità_ref e intensità odmr
    inizio, fine = trova_corrispondenza(intensità_ref, ch_A_od, mean=mean)
    # da questo indice parte, poi so che da lì ci sono altri 2799 punti
    
    # ritaglio gli spettri
    intensità_od = ch_A_od[inizio:fine]
    intensità_lock_in = ch_B_od[inizio:fine]
    
    return frequenze, intensità_od, intensità_lock_in, intensità_ref

def plotta_su_frequenze_odmr (reference: str, odmr: str, debug: bool = False, mean: bool = False):
    
    frequenze, intensità_od, intensità_lock_in, intensità_ref = leggi_files(reference, odmr, mean=True)
    
    # prima quella non lokkata
    plt.figure(figsize=(10, 3), dpi=200)
    plt.plot(frequenze, intensità_od, color='red')
    plt.xlabel('Frequenza [kHz]')
    plt.ylabel('Intensità PL [V]')
    plt.title('Intensità PL vs Frequenza')
    plt.grid(True)
    plt.show()
    
    if debug:
        # stampa le intensità di riferimento per vedere come si sovrappongono
        plt.figure(figsize=(10, 3), dpi=200)
        plt.plot(frequenze, intensità_ref, color='orange', label='REFERENCE')
        plt.xlabel('Frequenza [kHz]')
        plt.ylabel('Intensità PL [V]')
        plt.title('Intensità REFERENCE vs Frequenza')
        plt.grid(True)
        plt.show()
        
    # poi quella lokkata
    plt.figure(figsize=(10, 3), dpi=200)
    plt.plot(frequenze, intensità_lock_in, color='blue')
    plt.xlabel('Frequenza [kHz]')
    plt.ylabel('Intensità PL Lokkata [V]')
    plt.title('Intensità PL Lokkata vs Frequenza')
    plt.grid(True)
    plt.show()