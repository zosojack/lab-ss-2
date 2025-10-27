% filepath: prove_fit_6linee.m
% Fit di Exp.SampleFrame (e opzionalmente Field) per ottenere 6 deep in [2.5 3.2] GHz

clear; clc;

% --- EasySpin ---
espath = '/Users/zosojack/easyspin-6.0.11/easyspin';
if isempty(which('pepper')), addpath(genpath(espath)); rehash; end

% --- Sistema NV e Exp ---
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];        % [D E] MHz (E regola la separazione dei doppi)
Sys.lwpp = 6;             % restringi se serve separazione (es. 4–6 MHz)

Exp.Field = 16.85433;     % mT (verrà scalato leggermente opzionalmente)
Exp.mwRange = [2.5 3.2];  % FISSO
Exp.Harmonic = 0;

ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;  % fisso per NV in diamante
Exp.CrystalSymmetry = 227;

% --- Target (6 linee) ---
target = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];

% --- Angoli misurati B–NV (gradi) ---
theta_meas = [51 78 85];

% --- Assi NV cristallo (<111>) ---
Uall = [ 1  1  1;
         1 -1  1;
        -1  1  1;
        -1 -1  1] / sqrt(3);

% --- Stima iniziale di B (come prima) ---
subs = nchoosek(1:4,3); perms3 = perms(1:3);
best = struct('sse',Inf,'idx',[],'perm',[],'b',[]);
for si=1:size(subs,1)
    idx=subs(si,:); U=Uall(idx,:);
    for pi=1:size(perms3,1)
        c = cosd(theta_meas(perms3(pi,:))).';
        b = U\c; if norm(b)==0, continue; end
        b = b/norm(b);
        thp = acosd(U*b);
        sse = sum((thp - theta_meas(perms3(pi,:)).').^2);
        if sse<best.sse, best=struct('sse',sse,'idx',idx,'perm',perms3(pi,:),'b',b); end
    end
end

b0 = best.b(:);
phi0 = atan2(b0(2), b0(1));
theta0 = acos(b0(3));
Opt.Sites = best.idx;   % usa i 3 siti osservati

% --- Parametri estrazione deep ---
threshold = -0.1;   % ips negativo
dminGHz   = 0.030;  % distanza minima
promMin   = 0.06;   % prominenza minima
dx_fallbackGHz = 0.001; % stima passo se serve

get_deeps = @(x,ips) local_get_deeps(x,ips,threshold,dminGHz,promMin,dx_fallbackGHz);

% --- Piccolo fit locale: varia phi/theta (±5°) e, opzionale, Field (±2%) ---
dgrid = (-5:1:5)*pi/180;         % +/-5° a passi di 1°
fieldScales = 1;                 % tieni fisso il campo
% Se vuoi consentire piccola calibrazione, usa: fieldScales = 0.98:0.01:1.02;

bestFit = struct('sse',Inf,'phi',phi0,'theta',theta0,'fs',1,'deeps',[]);

for fs = fieldScales
    Exp.Field = 16.85433 * fs;
    for dph = dgrid
        for dth = dgrid
            Exp.SampleFrame = [phi0+dph, theta0+dth, 0];
            try
                [x,y] = pepper(Sys, Exp, Opt);
            catch
                continue;
            end
            ips = -y / max(abs(y)+eps);
            deeps = get_deeps(x,ips);
            if numel(deeps) < 6
                % penalizza se meno di 6: salta
                continue;
            end
            % match greedily 6 deeps ai 6 target
            sse = local_match_sse(deeps, target);
            if sse < bestFit.sse
                bestFit = struct('sse',sse,'phi',Exp.SampleFrame(1), ...
                    'theta',Exp.SampleFrame(2),'fs',fs,'deeps',deeps);
            end
        end
    end
end

% --- Simulazione finale con i migliori parametri trovati ---
Exp.Field = 16.85433 * bestFit.fs;
Exp.SampleFrame = [bestFit.phi, bestFit.theta, 0];
[x,y] = pepper(Sys, Exp, Opt);
ips = -y / max(abs(y)+eps);
deeps = get_deeps(x,ips);

fprintf('Best: phi=%.2f°, theta=%.2f°, FieldScale=%.3f, SSE=%.3g, #deeps=%d\n', ...
    bestFit.phi*180/pi, bestFit.theta*180/pi, bestFit.fs, bestFit.sse, numel(deeps));

% --- Plot ---
figure('Position',[100 100 1100 520]);
plot(x, ips, 'b-', 'LineWidth', 1.5); grid on; hold on;
xlabel('Frequenza (GHz)'); ylabel('Intensità (norm.)');
title(sprintf('Fit 6 linee | phi=%.1f°, theta=%.1f° | Field %.5f mT', ...
    bestFit.phi*180/pi, bestFit.theta*180/pi, 16.85433*bestFit.fs));
% deeps selezionati (solo i 6 abbinati ai target)
[~, matchIdx] = local_match_sse(deeps, target);
xline(deeps(matchIdx), 'r-', 'LineWidth', 0.8);
% target
for i=1:numel(target), xline(target(i),'k:'); end
legend('pepper','deep selezionati','target','Location','best');

% ===== Funzioni locali =====
function [deeps] = local_get_deeps(x,ips,threshold,dminGHz,promMin,dx_fallback)
    if numel(x)>1, dx = mean(diff(x)); else, dx = dx_fallback; end
    if exist('findpeaks','file')==2
        s = -ips;
        [~, locsX] = findpeaks(s, x, ...
            'MinPeakHeight', -threshold, ...
            'MinPeakDistance', dminGHz, ...
            'MinPeakProminence', promMin);
        deeps = sort(locsX(:).');
    else
        winGHz = 0.02;
        dminSamp = max(1, round(dminGHz/dx));
        winSamp  = max(1, round(winGHz/dx));
        cand = find(ips(2:end-1) < ips(1:end-2) & ips(2:end-1) < ips(3:end)) + 1;
        cand = cand(ips(cand) < threshold);
        keep = false(size(cand));
        for ii=1:numel(cand)
            i0=cand(ii); L=max(1,i0-winSamp); R=min(numel(ips),i0+winSamp);
            prom = max(ips(L:R)) - ips(i0);
            keep(ii) = prom >= promMin;
        end
        cand = cand(keep);
        [~,ord] = sort(ips(cand),'ascend'); cand = cand(ord);
        sel = [];
        for ii=1:numel(cand)
            if isempty(sel) || all(abs(x(cand(ii)) - x(sel)) >= dminGHz)
                sel(end+1)=cand(ii); %#ok<AGROW>
            end
        end
        deeps = sort(x(sel));
    end
end

function [sse, pickIdx] = local_match_sse(deeps, target)
    % Abbina greedily 6 deeps ai 6 target (1:1, per prossimità)
    deeps = deeps(:).'; target = target(:).';
    % se più di 6 deeps, usa i più vicini ai target
    pickIdx = zeros(1,numel(target));
    avail = true(1,numel(deeps));
    sse = 0;
    for k=1:numel(target)
        [~,idx] = min(abs(deeps - target(k)) + (~avail)*1e3);
        pickIdx(k) = idx; avail(idx)=false;
        sse = sse + (deeps(idx)-target(k))^2;
    end
end