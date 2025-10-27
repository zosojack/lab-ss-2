clear all; clc;

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

% ----- Range di angoli -----
step = 1;  % più piccolo = più preciso, ma più lento
alpha_range = 0:step:110;
beta_range  = 0:step:180;
gamma_range = 0:step:180;

% ----- File di output -----
outfile = 'angoli_validi_step1_from0_to110.txt';
fid = fopen(outfile, 'w');
fprintf(fid, 'alpha_deg\tbeta_deg\tgamma_deg\tsse\tsum_diff\n');

% ----- Setup progress -----
total_steps = length(alpha_range)*length(beta_range)*length(gamma_range);
current_step = 0;

% ----- Best tracking -----
best_sse = inf;
best_angles = [NaN, NaN, NaN];

% ----- Parametri selezione deep -----
threshold = -0.1;   % intensità minima (ips)
dminGHz   = 0.035;  % distanza minima tra deep in GHz
promMin   = 0.06;   % prominenza minima (unità di ips normalizzato)
winGHz    = 0.02;   % finestra per prominenza nel fallback manuale

for a = alpha_range
    for b = beta_range
        for g = gamma_range
            current_step = current_step + 1;

            % --- Aggiorna progresso ogni 1000 passi ---
            if mod(current_step, 1000) == 0
                progress = (current_step / total_steps) * 100;
                fprintf('\rProgresso: %6.2f%%', progress);
            end

            % --- Simulazione ---
            Exp.SampleFrame = [a b g]*pi/180;
            [x, y] = pepper(Sys, Exp, Opt);
            denom = max(abs(y));
            if denom==0 || ~isfinite(denom)
                ips = -y;
            else
                ips = -y / denom;
            end

            % --- Selezione deep: threshold + distanza minima + prominence ---
            locs = [];
            if exist('findpeaks','file')==2
                % Usa X per avere MinPeakDistance in GHz
                s = -ips; % minimi di ips -> massimi
                [~, locs_val] = findpeaks(s, x, ...
                    'MinPeakHeight', -threshold, ...
                    'MinPeakDistance', dminGHz, ...
                    'MinPeakProminence', promMin);
                locs = locs_val(:).';
            else
                % Fallback manuale
                dx = mean(diff(x));
                dminSamp = max(1, round(dminGHz/dx));
                winSamp  = max(1, round(winGHz/dx));

                cand = find(ips(2:end-1) < ips(1:end-2) & ips(2:end-1) < ips(3:end)) + 1;
                cand = cand(ips(cand) < threshold);
                if ~isempty(cand)
                    keep = false(size(cand));
                    for ii = 1:numel(cand)
                        i0 = cand(ii);
                        L = max(1, i0 - winSamp);
                        R = min(length(ips), i0 + winSamp);
                        prom = max(ips(L:R)) - ips(i0); % "profondità" locale
                        keep(ii) = prom >= promMin;
                    end
                    cand = cand(keep);

                    % Ordina per profondità (più negativi prima)
                    [~, ord] = sort(ips(cand), 'ascend');
                    cand = cand(ord);

                    % Selezione greedy con distanza minima in GHz
                    sel = [];
                    for ii = 1:numel(cand)
                        if isempty(sel) || all(abs(x(cand(ii)) - x(sel)) >= dminGHz)
                            sel(end+1) = cand(ii); %#ok<AGROW>
                        end
                    end
                    locs = x(sel);
                end
            end

            % --- Condizione rigorosa + SSE + somma differenze ---
            if numel(locs) == n_targets
                locs = sort(locs);
                diff_vector = locs(:) - target_peaks(:);
                sse = sum(diff_vector.^2);
                sum_diff = sum(diff_vector);

                % aggiorna best se migliore (scrive solo miglioramenti)
                if sse <= best_sse
                    fprintf(fid, '%.2f\t%.2f\t%.2f\t%.8g\t%.8g\n', a, b, g, sse, sum_diff);
                    best_sse = sse;
                    best_angles = [a, b, g];
                end
            end
        end
    en
end

fprintf('\n✅ Angoli salvati in "%s"\n', outfile);
fprintf('Migliori angoli: [%.2f, %.2f, %.2f] con SSE = %.6g\n', best_angles(1), best_angles(2), best_angles(3), best_sse);
fclose(fid);