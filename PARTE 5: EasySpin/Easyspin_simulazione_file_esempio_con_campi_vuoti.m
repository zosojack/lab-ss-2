
%%% 1. Scarica il software Easyspin dal link: https://easyspin.org/download.html
%%% 2. segui le istruzione per l'installazione

%%% Qui di seguito è riportato un esempio di FILE per simulare spettri di risonanza magnetica
%%% nella modalità "frequency sweep" per fissato valore e orientazione del campo magnetico esterno. 

%%%Per dettagli su come dichiarare le funzioni di Easyspin fare riferiemento
%%% al link: https://easyspin.org/easyspin/documentation/index.html


clear; clc;   clf;
% Importa dati sperimentali 

    % esempio:   fid = fopen('filename.txt');
    % etc....


% Definisci l'Hamiltoniana di SPIN          
Sys.S = ;  % 
Sys.g = ;  % assumre free electron g value: g=2.0
Sys.D = [  ];  % inserire i valori di zero-field splitting assiale, D e rombico, E in unità di misura di MHz

Sys.Nucs = ''; % definire gli isotopi nucleari (tabella degli isotopi: https://easyspin.org/easyspin/documentation/isotopetable.html) 
Sys.A = ; % interazione iperfine in MHz 

Sys.lwpp = ;  % definire la larghezza di riga in MHz


% Definire i parametri sperimentali 

Exp.Field = ;  % campo magnetico applicato in mT 
Exp.mwRange = [ ]; % definire intervallo di frequenze in GHz
Exp.Harmonic=; % 0 per spettro di assorbimento, 1 per la derivata prima

ma = 54.73561;  % magic angle, degrees  
Exp.MolFrame = [45 ma 0]*pi/180;  % definire l'orientazione (angoli di Eulero in radianti) dei centri NV rispetto al riferimento del cristallo  
Exp.CrystalSymmetry = '';  % diamond space group, for info and visuliazation visit: https://www.ccdc.cam.ac.uk/structures/?  
Exp.SampleFrame = [alfa beta gamma]*pi/180; % definire l'orientazione (angoli di Eulero) del cristallo rispetto al riferimento del laboratorio (cioè del campo magnetico esterno assumendo che B // z 


Opt.Sites =  []; % lasciare campo vuoto per includere tutti i centri NV non equivalenti 

% Calcola lo spettro
[x,y]=pepper(Sys,Exp,Opt); % calcola lo spettro e restituisce i valori nelle variabili x (frequenza) e y (intensita). 

% Grafica i risultati
plot(x,y);  % grafica lo spettro sperimentale e quello simulato. Se necessario aggiusta i parametri della simulazione per ottenere un buon accordo. 

figure;

% % Calcolo dei livelli energetici 
 ori = [ ]*pi/180;  % definire una orientazione del campo magnetico (https://easyspin.org/easyspin/documentation/levels.html)
 Be=0:0.1:5;   % definire un intervallo di campo magnetico 
 E = levels(Sys,ori,Be);
 plot(Be, E);
 shg;


% Opzionale: energy levels plot e transizioni 
% MF = 2.87;
% FieldRange = [0 500]; 
% levelsplot(Sys,ori,FieldRange,MF);
