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

% ----- Picchi target (MHz) -----
target_peaks = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];
n_targets = length(target_peaks);

% ----- Range di angoli -----
step = 1;  % più piccolo = più preciso, ma più lento
alpha_range = 172:step:180;
beta_range  = 0:step:180;
gamma_range = 0:step:180;

% ----- File di output -----
outfile = 'angoli_validi_step1_from172_to180.txt';
fid = fopen(outfile, 'w');
fprintf(fid, 'alpha_deg\tbeta_deg\tgamma_deg\tsse\tsum_diff\n');  % aggiunta colonna sum_diff


% ----- Setup progress -----
total_steps = length(alpha_range)*length(beta_range)*length(gamma_range);
current_step = 0;

% ----- Best tracking -----
best_sse = inf;
best_angles = [NaN, NaN, NaN];

for a = alpha_range
    for b = beta_range
        for g = gamma_range
            current_step = current_step + 1;

            % --- Aggiorna progresso ogni 360 passi ---
            if mod(current_step, 360) == 0
                progress = (current_step / total_steps) * 100;
                fprintf('\rProgresso: %6.2f%%', progress);
            end

            % --- Simulazione ---
            Exp.SampleFrame = [a b g]*pi/180;
            [x, y] = pepper(Sys, Exp, Opt);
            ips = -y / max(y);

            % --- Trova deep con threshold ---
            threshold = -0.1;
            locs = [];
            for i = 2:length(ips)-1
                if ips(i) < ips(i-1) && ips(i) < ips(i+1) && ips(i) < threshold
                    locs(end+1) = x(i);
                end
            end

            % --- Condizione rigorosa + SSE + somma differenze ---
            if numel(locs) == n_targets
                locs = sort(locs);  % garantisce stesso ordine dei target
                diff_vector = locs(:) - target_peaks(:);
                sse = sum(diff_vector.^2);  % scarto quadratico
                sum_diff = sum(diff_vector);  % somma semplice delle differenze
                
                % aggiorna best se migliore
                if sse <= best_sse
                    % salva angoli + sse + sum_diff con precisione aumentata
                    fprintf(fid, '%.2f\t%.2f\t%.2f\t%.8g\t%.8g\n', a, b, g, sse, sum_diff);
                    best_sse = sse;
                    best_angles = [a, b, g];
                end
            end
        end
    end
end

fprintf('\n✅ Angoli salvati in "%s"\n', outfile);
fprintf('Migliori angoli: [%.2f, %.2f, %.2f] con SSE = %.6g\n', best_angles(1), best_angles(2), best_angles(3), best_sse);
fclose(fid);
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

% ----- Picchi target (MHz) -----
target_peaks = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];
n_targets = length(target_peaks);

% ----- Range di angoli -----
step = 1;  % più piccolo = più preciso, ma più lento
alpha_range = 172:step:180;
beta_range  = 0:step:180;
gamma_range = 0:step:180;

% ----- File di output -----
outfile = 'angoli_validi_step1_from172_to180.txt';
fid = fopen(outfile, 'w');
fprintf(fid, 'alpha_deg\tbeta_deg\tgamma_deg\tsse\tsum_diff\n');  % aggiunta colonna sum_diff


% ----- Setup progress -----
total_steps = length(alpha_range)*length(beta_range)*length(gamma_range);
current_step = 0;

% ----- Best tracking -----
best_sse = inf;
best_angles = [NaN, NaN, NaN];

for a = alpha_range
    for b = beta_range
        for g = gamma_range
            current_step = current_step + 1;

            % --- Aggiorna progresso ogni 360 passi ---
            if mod(current_step, 360) == 0
                progress = (current_step / total_steps) * 100;
                fprintf('\rProgresso: %6.2f%%', progress);
            end

            % --- Simulazione ---
            Exp.SampleFrame = [a b g]*pi/180;
            [x, y] = pepper(Sys, Exp, Opt);
            ips = -y / max(y);

            % --- Trova deep con threshold ---
            threshold = -0.1;
            locs = [];
            for i = 2:length(ips)-1
                if ips(i) < ips(i-1) && ips(i) < ips(i+1) && ips(i) < threshold
                    locs(end+1) = x(i);
                end
            end

            % --- Condizione rigorosa + SSE + somma differenze ---
            if numel(locs) == n_targets
                locs = sort(locs);  % garantisce stesso ordine dei target
                diff_vector = locs(:) - target_peaks(:);
                sse = sum(diff_vector.^2);  % scarto quadratico
                sum_diff = sum(diff_vector);  % somma semplice delle differenze
                
                % aggiorna best se migliore
                if sse <= best_sse
                    % salva angoli + sse + sum_diff con precisione aumentata
                    fprintf(fid, '%.2f\t%.2f\t%.2f\t%.8g\t%.8g\n', a, b, g, sse, sum_diff);
                    best_sse = sse;
                    best_angles = [a, b, g];
                end
            end
        end
    end
end

fprintf('\n✅ Angoli salvati in "%s"\n', outfile);
fprintf('Migliori angoli: [%.2f, %.2f, %.2f] con SSE = %.6g\n', best_angles(1), best_angles(2), best_angles(3), best_sse);
fclose(fid);