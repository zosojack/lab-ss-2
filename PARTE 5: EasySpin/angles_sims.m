clear all; clc;

% ----- Sistema di spin -----
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

% ----- Parametri sperimentali -----
Exp.Field = 1.5;              % mT
Exp.mwRange = [2.5 3.2];      % GHz
Exp.Harmonic = 0;
ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

% ----- Target peaks (opzionale) -----
target_peaks = [2.816 2.844 2.874 2.899];

% ---------- Parametri tuning del peak finder ----------
threshold = -0.05;      % intensità minima richiesta
dminGHz   = 0.001;      % distanza minima tra picchi
promMin   = 0.02;       % prominenza minima
winGHz    = 0.02;       % finestra locale per prominenza (fallback)

% ======= INSERISCI QUI LA TERNA DA TESTARE =======
alpha = 90;
beta  = 45;
gamma = 30;
% ==================================================

fprintf("\n=== Test per terna (alpha=%.1f°, beta=%.1f°, gamma=%.1f°) ===\n", ...
        alpha, beta, gamma);

% ----- Simulazione -----
Exp.SampleFrame = [alpha beta gamma]*pi/180;
[x, y] = pepper(Sys, Exp, Opt);
ips = -y / max(abs(y));

% ----- Peak detection -----
locs = [];
promList = [];
depthList = [];

dx = mean(diff(x));
dminSamp = max(1, round(dminGHz/dx));
winSamp  = max(1, round(winGHz/dx));

% Trova minimi locali
cand = find(ips(2:end-1) < ips(1:end-2) & ips(2:end-1) < ips(3:end)) + 1;
cand = cand(ips(cand) < threshold);

for ii = 1:numel(cand)
    i0 = cand(ii);
    L = max(1, i0-winSamp);
    R = min(length(ips), i0+winSamp);
    prom = max(ips(L:R)) - ips(i0);   % prominenza grezza
    if prom >= promMin
        locs(end+1) = x(i0); %#ok<AGROW>
        promList(end+1) = prom;
        depthList(end+1) = ips(i0);
    end
end

% ----- Output diagnostico -----
fprintf("\nTrovati %d picchi.\n", numel(locs));
for k = 1:numel(locs)
    fprintf("  Picco %d: pos = %.4f GHz, depth = %.4f, prom = %.4f\n", ...
        k, locs(k), depthList(k), promList(k));
end
if isempty(locs)
    fprintf("  Nessun picco soddisfa i criteri.\n");
end

% ----- Plot diagnostico -----
figure; hold on; grid on;

plot(x, ips, 'LineWidth', 1.3);
yline(threshold, 'k--', 'LineWidth', 1.2);

% plotting target peaks
for p = target_peaks
    xline(p, 'r--', 'LineWidth', 1.2);
end

% plotting picchi trovati
if ~isempty(locs)
    plot(locs, depthList, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
end

xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title(sprintf('alpha=%.1f°, beta=%.1f°, gamma=%.1f°', alpha, beta, gamma));

legend('Spettro', 'Threshold', 'Target peaks', 'Picchi trovati');

hold off;