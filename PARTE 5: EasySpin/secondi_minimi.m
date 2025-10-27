% secondi_minimi.m
% Visualizza configurazioni di angoli (alpha,beta,gamma)
% Incolla uno o più array Nx3 in angleSets.

clear; clc;

% --- Assicura EasySpin nel path ---
espath = '/Users/zosojack/easyspin-6.0.11/easyspin';
if isempty(which('pepper'))
    if exist(espath,'dir')
        addpath(genpath(espath));
        rehash; rehash toolboxcache;
    end
end
if isempty(which('pepper'))
    fprintf('ERRORE: pepper non trovata.\nwhich -all pepper:\n');
    disp(which('-all','pepper'));
    error('Aggiungi EasySpin al path: %s', espath);
end

% --- Sistema di spin / Exp EasySpin ---
Sys.S = 1;
Sys.g = 2.0;
Sys.D = [2870 10];
Sys.Nucs = 'C';
Sys.A = [127, 127, 127];
Sys.lwpp = 8;

Exp.Field = 16.85433;       % mT
Exp.mwRange = [2.5 3.2];    % GHz
Exp.Harmonic = 0;
ma = 54.73561;
Exp.MolFrame = [45 ma 0]*pi/180;
Exp.CrystalSymmetry = 227;
Opt.Sites = [];

% --- Target (GHz) ---
target_peaks = [2.58303, 2.76875, 2.83447, 2.94180, 2.99815, 3.12104];
n_targets = numel(target_peaks);

% --- Parametri deep ---
threshold = -0.1;

% --- INCOLLA QUI I TUOI ARRAY Nx3 ---
anglesA = [
  9   60   16;
  9   60  106;
  9  120   74;
  9  120  164;
 10   60   16;
 10   60  106;
 10  120   74;
 10  120  164
];

anglesB = [
 110   60   74;
 110   60  164;
 110  120   16;
 110  120  106;
 111   60   74;
 111   60  164;
 111  120   16;
 111  120  106;
 112   60   74;
 112   60  164;
 112  120   16;
 112  120  106;
 113   60   74;
 113   60  164;
 113  120   16;
 113  120  106;
 114   60   74;
 114   60  164;
 114  120   16;
 114  120  106
];

angleSets = { anglesA; anglesB };  % aggiungi altri set se vuoi

if isempty(angleSets)
    error('Nessun set di angoli definito. Incolla uno o più array Nx3 in angleSets.');
end

% --- Una figura per ogni set, griglia con divisori di N ---
for s = 1:numel(angleSets)
    angles = angleSets{s};
    if size(angles,2) < 3
        warning('Set %d ignorato: non ha 3 colonne.', s);
        continue;
    end
    nShow = size(angles,1);

    % Griglia: fattorizza N in righe x colonne (la più "quadrata" possibile)
    r = floor(sqrt(nShow));
    while r > 1 && mod(nShow,r) ~= 0
        r = r - 1;
    end
    if mod(nShow,r)==0
        nRows = r; nCols = nShow/r;
    else
        nRows = 1; nCols = nShow;   % N primo
    end

    % Precalcolo spettri e limiti
    spectra = cell(nShow,1);
    xs = cell(nShow,1);
    xlimits = [];
    ylimits = [Inf, -Inf];

    for k = 1:nShow
        a = angles(k,1); b = angles(k,2); g = angles(k,3);
        Exp.SampleFrame = [a b g]*pi/180;
        [x,y] = pepper(Sys, Exp, Opt);

        denom = max(abs(y));
        if ~isfinite(denom) || denom==0
            ips = -y;
        else
            ips = -y/denom;
        end

        xs{k} = x;
        spectra{k} = ips;

        if isempty(xlimits)
            xlimits = [min(x), max(x)];
        end
        ylimits(1) = min(ylimits(1), min(ips));
        ylimits(2) = max(ylimits(2), max(ips));
    end

    % Figura per il set corrente
    figure('Name',sprintf('Set %d (%d config) - %dx%d',s,nShow,nRows,nCols), ...
           'NumberTitle','off','Position',[100 100 1400 900]);
    tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');

    for j = 1:nShow
        nexttile;
        plot(xs{j}, spectra{j}, 'b-','LineWidth',1); hold on;

        % Linee verticali target
        for tt = 1:n_targets
            xline(target_peaks(tt),'g--','LineWidth',0.8);
        end

        % Deep + metriche
        ips = spectra{j}; x = xs{j};
        locs = [];
        for i = 2:length(ips)-1
            if ips(i) < ips(i-1) && ips(i) < ips(i+1) && ips(i) < threshold
                locs(end+1) = x(i); %#ok<SAGROW>
            end
        end
        sse = NaN; sumdiff = NaN;
        if numel(locs) == n_targets
            locs = sort(locs);
            diffv = locs(:) - target_peaks(:);
            sse = sum(diffv.^2);
            sumdiff = sum(diffv);
        end

        title(sprintf('(%d/%d) a=%.0f b=%.0f g=%.0f | SSE=%s | SUMDIFF=%s', ...
              j, nShow, angles(j,1), angles(j,2), angles(j,3), ...
              num2str(sse,'%.6g'), num2str(sumdiff,'%.6g')), 'FontSize',9);

        xlim(xlimits);
        if all(isfinite(ylimits)) && ylimits(1) < ylimits(2)
            ypad = 0.05*(ylimits(2)-ylimits(1));
            ylim([ylimits(1)-ypad, ylimits(2)+0.02*(ylimits(2)-ylimits(1))]);
        else
            axis tight;
        end

        xlabel('GHz'); ylabel('int norm'); grid on;
    end
end

disp('Plot completati.');