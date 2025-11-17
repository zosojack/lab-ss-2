clear all; clc;

% ----- Sistema di spin -----
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

% ----- Parametri sperimentali -----
Exp.Field = 1.5;       % mT
Exp.mwRange = [2.5 3.2];    % GHz
Exp.Harmonic = 0;
ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

% ----- Target peaks -----
target_peaks = [2.816 2.844 2.874 2.899];

% ======= QUI IMPOSTI I TUOI ANGOLI =======
alpha = 90;    % <-- N1
beta  = 45;    % <-- N2

gamma_i   = 10;   % <-- N3_i
gamma_step = 5;
N = 10;           % numero di grafici
% gamma_f = gamma_i + gamma_step * N
% ==========================================

threshold = -0.1;

for k = 0:N-1
    gamma = gamma_i + k * gamma_step;

    % ----- Simulazione -----
    Exp.SampleFrame = [alpha beta gamma] * pi/180;
    [x, y] = pepper(Sys, Exp, Opt);
    ips = -y / max(y);

    % ----- Trova picchi in base al tuo criterio -----
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

    % threshold line
    yline(threshold, 'k--', 'LineWidth', 1.2);

    % target peaks
    for p = target_peaks
        xline(p, 'r--', 'LineWidth', 1.2);
    end

    % punti trovati dal detector
    if ~isempty(locs)
        plot(locs, interp1(x, ips, locs), 'bo', 'MarkerSize', 7, 'LineWidth', 1.5);
    end

    xlabel('Frequenza (GHz)');
    ylabel('Intensità normalizzata');
    title(sprintf('alpha = %.1f°, beta = %.1f°, gamma = %.1f°', alpha, beta, gamma));
    legend('Spettro', 'Threshold', 'Target peaks', 'Picchi simulati');

    hold off;

    % stampa info a terminale
    fprintf('gamma = %.1f ---> picchi trovati: ', gamma);
    disp(locs);
end