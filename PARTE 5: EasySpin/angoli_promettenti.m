% Visualizza configurazioni di angoli (alpha,beta,gamma)
% Legge angoli da angoli_promettenti.mat (variabile 'angles' Nx3) o .txt.

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

% --- Lettura angoli ---
angles = [];
try
    if exist('angoli_promettenti.mat','file')
        tmp = load('angoli_promettenti.mat');
        if isfield(tmp,'angles') && size(tmp.angles,2)>=3
            angles = tmp.angles(:,1:3);
        end
    end
catch
    angles = [];
end
if isempty(angles) && exist('angoli_promettenti.txt','file')
    angles = readmatrix('angoli_promettenti.txt');
    angles = angles(:,1:3);
end
if isempty(angles)
    error('Nessun angolo trovato. Fornisci angoli in angoli_promettenti.mat (variabile angles Nx3) o .txt.');
end

% --- Configurazione visualizzazione ---
nShow = size(angles,1);
nPerFig = 36;                 % usa 9 per 3x3 per figura
if nPerFig==36, nRows=6; nCols=6;
elseif nPerFig==9, nRows=3; nCols=3;
else, nRows = ceil(sqrt(nPerFig)); nCols = ceil(nPerFig/nRows);
end

% --- Parametri di plotting ---
threshold = -0.1;            % deep threshold
xlimits = [];                % verrà impostato dal primo spettro
ylimits = [Inf, -Inf];       % [min, max] globale

% --- Precalcolo spettri ---
spectra = cell(nShow,1);
xs = cell(nShow,1);
for k = 1:nShow
    a = angles(k,1); b = angles(k,2); g = angles(k,3);
    Exp.SampleFrame = [a b g]*pi/180;
    [x,y] = pepper(Sys, Exp, Opt);

    % Normalizzazione robusta (evita divisione per 0)
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

% --- Plot a pagine ---
nPages = ceil(nShow / nPerFig);
for p = 1:nPages
    fig = figure('Name',sprintf('Configurazioni %d/%d',p,nPages),'NumberTitle','off',...
                 'Position',[100 100 1400 900]);
    tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');
    startIdx = (p-1)*nPerFig + 1;
    endIdx = min(p*nPerFig, nShow);

    for j = startIdx:endIdx
        nexttile;
        plot(xs{j}, spectra{j}, 'b-','LineWidth',1); hold on;

        % Linee verticali target
        for tt = 1:n_targets
            xline(target_peaks(tt),'g--','LineWidth',0.8);
        end

        % Deep locali + metriche
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

        title(sprintf('(%d) a=%.0f b=%.0f g=%.0f\nSSE=%s SUMDIFF=%s', ...
              j, angles(j,1), angles(j,2), angles(j,3), ...
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

    % Riempi tile vuoti se la pagina non è piena
    nTilesUsed = endIdx - startIdx + 1;
    for tfill = nTilesUsed+1:nPerFig
        nexttile; axis off;
    end

    if p < nPages
        resp = input('Invio per pagina successiva, q per uscire: ','s');
        if strcmpi(strtrim(resp),'q')
            break;
        end
        close(fig);
    end
end

disp('Plot completati.');
