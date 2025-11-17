clear all; clc;

%% ============================================================
%% Parametri del sistema (come il tuo esempio)
%% ============================================================
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

Exp.mwRange = [2.5 3.2];   % GHz
Exp.Harmonic = 0;

ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

%% ============================================================
%% Angoli da usare (fissi)
%% ============================================================
alpha = 10;
beta  = 66;
gamma = 49;

Exp.SampleFrame = [alpha beta gamma] * pi/180;

%% ============================================================
%% Rampa dei campi (mT)
%% ============================================================
campi = [ ...
    0.00147035
    0.00255966
    0.00345788
    0.00671791
    0.00857636
    0.01110222
    0.01270503
    0.01685433];

%% ============================================================
%% MATRICE DEI PICCHI (DA COMPILARE)
%% Ogni riga = picchi sperimentali per quel campo
%% Usa NaN per celle vuote
%% ============================================================
picchi_matrix = {
    [2.760 2.816 2.844 2.874 2.899 2.969]    % campo 1 → qui metti i tuoi picchi
    [2.738 2.793 2.833 2.841 2.872 2.881 2.917 2.986]    % campo 2
    [2.723 2.778 2.826 2.829 2.837 2.879 2.888 2.892 2.933 3.002]    % ...
    [2.669 2.724 2.792 2.807 2.832 2.895 2.919 2.934 2.987 3.055]
    [2.643 2.698 2.766 2.810 2.826 2.898 2.912 2.947 3.003 3.071]
    [2.604 2.660 2.729 2.789 2.830 2.925 2.961 3.001 3.056 3.126]
    [2.587 2.643 2.713 2.790 2.840 2.925 2.968 3.019 3.074 3.145]
    [ ]
};

%% ============================================================
%% Loop su tutti i campi
%% ============================================================
threshold = 0.1;   % per trovare deep simulati

for k = 1:length(campi)
    
    fprintf("\n=== Campo = %.6f mT ===\n", campi(k));
    
    Exp.Field = campi(k)*1000;   % EasySpin vuole mT → moltiplico per 1000
    
    % --- Simulazione ---
    [x, y] = pepper(Sys, Exp, Opt);
    ips = -y / max(y);
    
    % --- Trova picchi simulati ---
    locs = [];
    for i = 2:length(ips)-1
        if ips(i) < ips(i-1) && ips(i) < ips(i+1) && ips(i) < threshold
            locs(end+1) = x(i);
        end
    end
    locs = sort(locs);
    
    % --- Recupera i picchi sperimentali ---
    picchi_exp = picchi_matrix{k};
    
    %% --- Plot ---
    figure('Name',sprintf('Campo %.6f mT',campi(k))); hold on; grid on;
    
    plot(x, ips, 'LineWidth', 1.2);  % spettro
    
    % --- linee verticali tratteggiate per i picchi sperimentali ---
    for p = picchi_exp
        xline(p, 'r--', 'LineWidth', 1.2);
    end
    
    % --- picchi simulati ---
    plot(locs, interp1(x,ips,locs), 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    xlabel('Frequenza (GHz)');
    ylabel('Intensità normalizzata');
    title(sprintf('Campo = %.6f mT  |  Frame = [%d %d %d]°', ...
          campi(k), alpha, beta, gamma));
    
    legend('Simulato','Picchi sperimentali','Picchi trovati');
    
    fprintf("Picchi simulati trovati:\n");
    disp(locs(:));
end