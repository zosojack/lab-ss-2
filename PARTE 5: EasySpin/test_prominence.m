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

% ----- Picchi target (GHz) -----
target_peaks = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];
n_targets = length(target_peaks);

% ----- Angoli da testare -----
a = 110;  % alpha in gradi
b = 33;   % beta in gradi
g = 153;  % gamma in gradi

% ----- Simulazione -----
Exp.SampleFrame = [a b g]*pi/180;
[x, y] = pepper(Sys, Exp, Opt);
ips = -y / max(abs(y));  % normalizzazione robusta

% ----- Trova tutti i deep (minimi) locali (senza filtri) -----
all_deeps_x = [];
all_deeps_y = [];
for i = 2:length(ips)-1
    if ips(i) < ips(i-1) && ips(i) < ips(i+1)
        all_deeps_x(end+1) = x(i);
        all_deeps_y(end+1) = ips(i);
    end
end

% ----- Parametri filtri deep -----
threshold = -0.1;   % intensità minima (ips)
dminGHz   = 0.035;  % distanza minima tra deep in GHz
promMin   = 0.06;   % prominenza minima (in unità di ips)
winGHz    = 0.02;   % finestra per prominenza nel fallback manuale

dx = mean(diff(x));
dminSamp = max(1, round(dminGHz/dx));
winSamp  = max(1, round(winGHz/dx));

% ----- Selezione deep con threshold + distanza minima + prominence -----
sel_x = []; sel_y = [];
if exist('findpeaks','file')==2
    % Usa findpeaks su -ips (i minimi di ips diventano massimi)
    s = -ips;
    [~, locs_idx] = findpeaks(s, ...
        'MinPeakHeight', -threshold, ...      % ips < threshold  -> s > -threshold
        'MinPeakDistance', dminSamp, ...      % distanza minima in campioni
        'MinPeakProminence', promMin);        % prominenza minima
    sel_x = x(locs_idx);
    sel_y = ips(locs_idx);
else
    % Fallback manuale: minimi locali + prominenza locale + distanza minima
    cand = find(ips(2:end-1) < ips(1:end-2) & ips(2:end-1) < ips(3:end)) + 1;
    cand = cand(ips(cand) < threshold);
    if ~isempty(cand)
        keep = false(size(cand));
        for ii = 1:numel(cand)
            i0 = cand(ii);
            L = max(1, i0 - winSamp);
            R = min(length(ips), i0 + winSamp);
            prom = max(ips(L:R)) - ips(i0);   % "profondità" locale
            keep(ii) = prom >= promMin;
        end
        cand = cand(keep);

        % Ordina per profondità (più negativi prima)
        [~, ord] = sort(ips(cand), 'ascend');
        cand = cand(ord);

        % Selezione greedy con distanza minima
        sel_idx = [];
        for ii = 1:numel(cand)
            if isempty(sel_idx) || all(abs(x(cand(ii)) - x(sel_idx)) >= dminGHz)
                sel_idx(end+1) = cand(ii); %#ok<AGROW>
            end
        end
        sel_x = x(sel_idx);
        sel_y = ips(sel_idx);
    end
end

% ----- Visualizza lo spettro e i deep -----
figure('Position', [100, 100, 1200, 800]);

% Grafico principale dello spettro
subplot(2,1,1);
plot(x, ips, 'b-', 'LineWidth', 1.5); hold on;
% Minimi locali (tutti)
scatter(all_deeps_x, all_deeps_y, 36, 'ro');
% Minimi selezionati (filtrati)
scatter(sel_x, sel_y, 60, 'g^', 'filled');
xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title(sprintf('Spettro simulato [\\alpha=%.1f°, \\beta=%.1f°, \\gamma=%.1f°]', a, b, g));
grid on;
for i = 1:length(target_peaks)
    xline(target_peaks(i), 'g--', 'LineWidth', 1);
end
legend('Spettro', 'Minimi locali', 'Selezionati (thr+dist+prom)', 'Target', 'Location', 'best');

% Grafico con vari threshold (riapplica anche distanza e prominence)
subplot(2,1,2);
plot(x, ips, 'b-', 'LineWidth', 1.2); hold on;
thresholds = [-0.05, -0.1, -0.15, -0.2, -0.3, -0.4];
colors = lines(length(thresholds));

for t = 1:length(thresholds)
    th = thresholds(t);
    deeps_x = []; deeps_y = [];
    if exist('findpeaks','file')==2
        s = -ips;
        [~, idx] = findpeaks(s, ...
            'MinPeakHeight', -th, ...
            'MinPeakDistance', dminSamp, ...
            'MinPeakProminence', promMin);
        deeps_x = x(idx);
        deeps_y = ips(idx);
    else
        cand = find(ips(2:end-1) < ips(1:end-2) & ips(2:end-1) < ips(3:end)) + 1;
        cand = cand(ips(cand) < th);
        if ~isempty(cand)
            keep = false(size(cand));
            for ii = 1:numel(cand)
                i0 = cand(ii);
                L = max(1, i0 - winSamp);
                R = min(length(ips), i0 + winSamp);
                prom = max(ips(L:R)) - ips(i0);
                keep(ii) = prom >= promMin;
            end
            cand = cand(keep);
            [~, ord] = sort(ips(cand), 'ascend');
            cand = cand(ord);
            sel_idx = [];
            for ii = 1:numel(cand)
                if isempty(sel_idx) || all(abs(x(cand(ii)) - x(sel_idx)) >= dminGHz)
                    sel_idx(end+1) = cand(ii); %#ok<AGROW>
                end
            end
            deeps_x = x(sel_idx);
            deeps_y = ips(sel_idx);
        end
    end
    scatter(deeps_x, deeps_y, 48, colors(t,:), 'filled');
    % Annotazione del conteggio
    text(3.15, th, sprintf('thr=%.2f: %d deep', th, numel(deeps_x)), 'Color', colors(t,:), 'FontSize', 9);
    yline(th, '--', 'Color', colors(t,:), 'LineWidth', 1);
end

xlabel('Frequenza (GHz)');
ylabel('Intensità normalizzata');
title('Selezione deep con soglia, distanza minima e prominenza');
grid on;

% Stampa informazioni
fprintf('Minimi locali trovati: %d\n', numel(all_deeps_x));
fprintf('Minimi selezionati (thr=%.2f, dmin=%.3f GHz, prom>=%.3f): %d\n', threshold, dminGHz, promMin, numel(sel_x));
if numel(sel_x)==n_targets
    sel_sorted = sort(sel_x);
    diffv = sel_sorted(:) - target_peaks(:);
    sse = sum(diffv.^2);
    sumdiff = sum(diffv);
    fprintf('SSE=%.6g, SUMDIFF=%.6g\n', sse, sumdiff);
end