clear all; clc;

% ----- Sistema di spin -----
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

% ----- Parametri sperimentali -----
Exp.Field = 1.47;            % mT
Exp.mwRange = [2.5 3.2];     % GHz
Exp.Harmonic = 0;
ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

% ----- Target peaks -----
target_peaks = [2.76 2.816 2.844 2.874 2.899 2.969];

% ============================================
%     QUI INSERISCI LA TERNA DA TESTARE
% ============================================
alpha = 11;   % <-- Sostituisci
beta  = 62;   % <-- Sostituisci
gamma = 49;   % <-- Sostituisci
% ============================================

Exp.SampleFrame = [alpha beta gamma]*pi/180;

% ----- Simulazione -----
[x, y] = pepper(Sys, Exp, Opt);
ips = -y / max(y);   % normalizzazione coerente con il tuo script

% ----- Trova deep sotto threshold = 0 -----
threshold = 0.1;
locs = [];
for i = 2:length(ips)-1
    if ips(i) < ips(i-1) && ips(i) < ips(i+1) && ips(i) < threshold
        locs(end+1) = x(i);
    end
end
locs = sort(locs);

% ----- Plot -----
figure; hold on; grid on;
plot(x, ips, 'LineWidth', 1.2);

% target peaks (linee rosse)
for p = target_peaks
    xline(p, 'r--', 'LineWidth', 1.2);
end

% picchi trovati (pallini blu)
plot(locs, interp1(x, ips, locs), 'bo', 'MarkerSize', 6, 'LineWidth', 1.5);

xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title(sprintf('alpha = %.1f°, beta = %.1f°, gamma = %.1f°', alpha, beta, gamma));
legend('Spettro simulato', 'Target peaks', 'Picchi trovati');

hold off;

% stampa informazioni a schermo
disp('Picchi simulati trovati:');
disp(locs(:));