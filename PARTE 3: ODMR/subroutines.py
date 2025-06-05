import numpy as np
import matplotlib.pyplot as plt

# rimuove i valori troppo distanti che non hanno nessun significato
def pulisci_valori_nosense (intensità):
    # prendo l'asse della frequenza, lo divido in quinti, prendo la regione centrale
    quinto = len(intensità) // 5
    intensita_centrale = intensità[2*quinto:3*quinto]
    
    # Calcolo media e deviazione standard della regione centrale
    mu = np.mean(intensita_centrale)
    sigma = np.std(intensita_centrale)
    
    # Copio l'intensità per poterla modificare
    intensita_pulita = np.copy(intensità)
    
    # Sostituisco i valori che distano più di 3 sigma dalla media con il valore precedente
    # lo faccio dal centro all'esterno
    for i in range(round(len(intensita_pulita)/2), len(intensita_pulita)):
        if abs(intensita_pulita[i] - intensita_pulita[i-1]) > 4 * sigma:
            intensita_pulita[i] = intensita_pulita[i-1]
            
    # devo fare la stessa cosa anche nell'altro verso
    for i in range(int(len(intensita_pulita)/2), -1, -1):
        if abs(intensita_pulita[i+1] - intensita_pulita[i]) > 4 * sigma:
            intensita_pulita[i] = intensita_pulita[i+1]    
    
    return intensita_pulita


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


from scipy.signal import find_peaks

def trova_valli(segnale, prominence=0.003, distance=10, width=0.1):
    """
    Trova le valli in un segnale usando find_peaks.
    
    Args:
        segnale: array del segnale da analizzare
        prominence: prominenza minima delle valli
        distance: distanza minima tra valli consecutive
        width: larghezza minima delle valli
    
    Returns:
        indici: posizioni delle valli nel segnale originale
        proprietà: caratteristiche delle valli trovate
    """
    # Inverto il segnale per trasformare le valli in picchi
    segnale_invertito = -segnale
    
    # Trovo i picchi nel segnale invertito (valli nell'originale)
    indici_valli, proprietà_valli = find_peaks(
        segnale_invertito,
        prominence=prominence,  # Quanto deve essere profonda la valle
        distance=distance,      # Distanza minima tra valli consecutive
        width=width             # Larghezza minima della valle
    )
    
    return indici_valli, proprietà_valli

# zero inteso come punto in cui la der prima inverte segno,
# siccome è traslata spesso non è proprio 0 ma -0.0005
def i_primo_zero (vettore, partenza, zero=0):
    first = partenza
    taglio = vettore[partenza:]
    for j in range(len(taglio)-2):
        if taglio[j] < zero and taglio[j+1] > zero or taglio[j] > zero and taglio[j+1] < zero:
            first = j + partenza  # restituisce l'indice del primo zero trovato
            return first
    return first  # se non trova nessun zero, ritorna l'indice di partenza

def trova_primi_zeri (spettro, indici_valli):
    # il segnale del lock in è una derivata prima! quindi la valle in realtà è
    # il primo zero che segue la rispettiva valle nel grafico lock-in!
    # non è precisamente zero però, è leggermente traslata sulle y
    # devo usare una media 
    
    indici_primi_zeri = []
    zero = np.nanmean(np.array(spettro))
    
    print("zero:", zero)
    
    # per ogni valle
    for i in indici_valli:
        # cerco il primo zero successivo
        indici_primi_zeri.append(i_primo_zero(spettro, partenza=i, zero=zero))
            
    return indici_primi_zeri


def migliora_corrispondenza (ref, lock_in):
    # devo pulire completamente il ref lasciando solo i picchi
    # quindi fare un confronto tra i picchi del ref e quelli del lock-in
    # quest'ultimo è comodo perché ha già i picchi ben definiti
    
    # valli ref
    i_valli_ref, _ = trova_valli(ref)
    i_valli_lock_in, _ = trova_valli(lock_in, prominence=0.00001, distance=10, width=0.0001)
    i_valli_od = trova_primi_zeri(lock_in, i_valli_lock_in)
    
    print("VALLI INDIVIDUATE:")
    print('ref', i_valli_ref)
    print('lock', i_valli_lock_in)
    print('od', i_valli_od)
    
    
    # cerca quelli più vicini
    min_distanza = 1000000
    indice_ref_minimo = -1
    indice_od_minimo = -1
    
    # METODO 1: per ogni valle nel riferimento
    '''
    for idx_ref in i_valli_ref:
        # Calcola la distanza da ogni valle nel lock-in
        for idx_od in i_valli_od:
            distanza = abs(idx_ref - idx_od)
            
            if distanza < min_distanza:
                min_distanza = distanza
                indice_ref_minimo = idx_ref
                indice_od_minimo = idx_od
    '''     
    # METODO 2: PER LA MINIMA DELLE VALLI (così sei sicuro che lo sia)
    # indice (nel vettore degli indici minimi) dell'indice che dà la valle minima
    indici_ordinati = np.argsort(ref[i_valli_ref])
    # indice che dà la valle minima:
    if i_valli_ref[indici_ordinati[0]] > 100:
        idx_ref = i_valli_ref[indici_ordinati[0]]
    elif i_valli_ref[indici_ordinati[1]] > 100:
        idx_ref = i_valli_ref[indici_ordinati[1]]
    else: 
        idx_ref = i_valli_ref[indici_ordinati[2]]
        
        
    print("minimo:", ref[idx_ref], "| indice:", idx_ref)
    
    # Calcola la distanza tra questa e ogni valle nel lock-in
    for idx_od in i_valli_od:
        distanza = abs(idx_ref - idx_od)
        
        if distanza < min_distanza:
            min_distanza = distanza
            indice_ref_minimo = idx_ref
            indice_od_minimo = idx_od
                
    # stimo l'offset con questa distanza!
    offset = (indice_od_minimo - indice_ref_minimo)
    
    return offset
    

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
    
    
def leggi_files (reference: str, odmr: str, index, mean):
    
    # leggo il reference
    data = np.loadtxt(reference, skiprows=9, delimiter=',')
    time_ref = data[:, 0]  # in s
    ch_A_ref = data[:, 1]  # in V - PL intensity
    ch_B_ref = data[:, 2]  # in V - reference

    # trovo la regione dello sweep e la ritaglio
    frequenze, intensità_ref = intensità_e_frequenze_sweep(time_ref, ch_A_ref, ch_B_ref)
    
    # pulisco i valori che non hanno senso
    intensità_ref_pulita = pulisci_valori_nosense(intensità_ref)
    
    # leggo odmr
    data2 = np.loadtxt(odmr, skiprows=9, delimiter=',')
    ch_A_od = data2[:, 1]  # in V - PL intensity
    ch_B_od = data2[:, 2]  # in V - lock-in
    
    # pulisco i valori che non hanno senso
    ch_A_od_pulita = pulisci_valori_nosense(ch_A_od)
    
    # cerco corrispondenzza tra intensità_ref e intensità odmr
    inizio, fine = trova_corrispondenza(intensità_ref_pulita, ch_A_od_pulita, mean=mean)
    
    # shifto per migliorare
    offset = migliora_corrispondenza(intensità_ref_pulita, ch_B_od[inizio:fine])
    
    if index == 0:
        offset -= 16
        
    elif index == 1:
        offset -= 46
        
    elif index == 4:
        offset -= 25
        
    elif index == 5:
        offset -= 2
        
    elif index == 6:
        offset -= 13
        
    elif index == 7:
        offset -= 35
        
    elif index == 8:
        offset -= 12

    inizio_shifted = inizio + offset
    fine_shifted = fine + offset
    
    print(f"INIZIO: {inizio}")
    print(f"FINE: {fine}")
    print(f"OFFSET: {offset}")
    
            
    # ritaglio gli spettri
    intensità_od = ch_A_od_pulita[inizio_shifted:fine_shifted]
    intensità_lock_in = ch_B_od[inizio_shifted:fine_shifted]
    
    return frequenze, intensità_od, intensità_lock_in, intensità_ref_pulita

def plotta_su_frequenze_odmr (reference: str, odmr: str, debug: bool = False, mean: bool = False, index=0):
    
    frequenze, intensità_od, intensità_lock_in, intensità_ref = leggi_files(reference, odmr, index=index, mean=True, )
    
    frequenze = np.array(frequenze)
    intensità_od = np.array(intensità_od)
    intensità_lock_in = np.array(intensità_lock_in)
    intensità_ref = np.array(intensità_ref)
    
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
        
        intensità_ref = np.array(intensità_ref)
        indici_valli, _ = trova_valli(intensità_ref)
        if len(indici_valli) == 0:
            print("Nessuna valle trovata nel segnale.")
        else:
            plt.scatter(frequenze[indici_valli], intensità_ref[indici_valli], color='green', label='Valli trovate')

        
        
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