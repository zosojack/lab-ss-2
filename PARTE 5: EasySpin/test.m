clear all; clc;

% ----- Aggiungi EasySpin al path -----
addpath('/Users/zosojack/easyspin-6.0.11/easyspin');

% ----- Sistema di spin -----
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

% ----- Parametri sperimentali -----
Exp.Field = 16.85433;       % mT
Exp.mwRange = [2.5 3.2];    % GHz
Exp.Harmonic = 0;
ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

% ----- Picchi target (MHz) -----
target_peaks = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];
n_targets = length(target_peaks);

% ----- Angoli da testare -----
a = 110;  % alpha in gradi
b = 66;  % beta in gradi
g = 169;  % gamma in gradi

% ----- Simulazione -----
Exp.SampleFrame = [a b g]*pi/180;
[x, y] = pepper(Sys, Exp, Opt);
ips = -y / max(y);

% ----- Trova tutti i deep (minimi) locali -----
all_deeps_x = [];
all_deeps_y = [];
for i = 2:length(ips)-1
    if ips(i) < ips(i-1) && ips(i) < ips(i+1)  % Cambiato > a <
        all_deeps_x(end+1) = x(i);
        all_deeps_y(end+1) = ips(i);
    end
end

% ----- Visualizza lo spettro e i deep -----
figure('Position', [100, 100, 1200, 800]);

% Grafico principale dello spettro
subplot(2,1,1);
plot(x, ips, 'b-', 'LineWidth', 1.5);
hold on;
scatter(all_deeps_x, all_deeps_y, 50, 'ro', 'filled');
xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title(sprintf('Spettro simulato con angoli [α=%.1f°, β=%.1f°, γ=%.1f°]', a, b, g));
grid on;

% Aggiungi linee verticali per i target
for i = 1:length(target_peaks)
    xline(target_peaks(i), 'g--', 'LineWidth', 1);
end
legend('Spettro', 'Deep (minimi) locali', 'Target', 'Location', 'best');

% Grafico con vari threshold
subplot(2,1,2);
plot(x, ips, 'b-', 'LineWidth', 1.5);
hold on;

% Testa diverse soglie negative
thresholds = [-0.05, -0.1, -0.15, -0.2, -0.3, -0.4];
colors = lines(length(thresholds));

for t = 1:length(thresholds)
    th = thresholds(t);
    % Trova deep che sono minori della soglia (più negativi)
    deeps_x = [];
    deeps_y = [];
    for i = 1:length(all_deeps_x)
        if all_deeps_y(i) < th
            deeps_x(end+1) = all_deeps_x(i);
            deeps_y(end+1) = all_deeps_y(i);
        end
    end
    
    % Plotta i deep sotto la soglia
    scatter(deeps_x, deeps_y, 50, colors(t,:), 'filled');
    
    % Stampa quanti deep sono stati trovati con questa soglia
    text(3.15, th, sprintf('threshold=%.2f: %d deep', th, length(deeps_x)), 'Color', colors(t,:));
    
    % Linea orizzontale per la soglia
    yline(th, '--', 'Color', colors(t,:), 'LineWidth', 1);
end

xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title('Confronto tra diverse soglie per il rilevamento dei deep');
grid on;

% Stampa informazioni sui deep e valori
disp('Tutti i deep locali trovati:');
for i = 1:length(all_deeps_x)
    fprintf('Deep a %.5f GHz, intensità: %.5f\n', all_deeps_x(i), all_deeps_y(i));
end