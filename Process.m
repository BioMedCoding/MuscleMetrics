%% ========================== Note iniziali ================================
% Allo stato attuale la stima della CV risulta essere errata, si ha un
% canale in particolare che riprota valori non fisiologici



%% ========================== Inizializzazione =============================
clear
close all
clc



%% ========================= Parametri generali ============================
mostra_plot_import = true;
mostra_plot_differenziali = true;

mostra_fatigue_normalizzati_singoli = true;
mostra_fatigue_normalizzati__medi = true;
mostra_fatigue_singoli = true;
mostra_fatigue_medi = true;


separa_muscoli = true;      % Se vero, va a passare alla funzione di analisi solamente la porzione di analisi definita dalla variabile is_bicipite
is_bicipite = true;         % Se vero, si passano all'analisi solamente i primi 8 canali della matrice, altrimenti gli ultimi 8
canale_plot = 1;            % Fino a risoluzione, perche non dia problemi deve essere tra 1 e 4
offset_canali = 8;          % Colonna che separa primo e secondo gruppo muscolare


risoluzione_calcolo = 0.5;  % Risoluzione del fatigue plot, in secondi


nome_segnale_calcolo_metriche = 'sig_singolo_diff';     % Nome del segnale che si vuole utilizzare nella parte dedicata ai fatigue plot per le metriche di ampiezza e frequenza
                                % Le possibilità sono: 'sig_mono', 'sig_singolo_diff', 'sig_doppio_diff'

nome_segnale_calcolo_cv = 'sig_singolo_diff';   % Nome del segnale che si vuole utilizzare nella parte dedicata ai fatigue plot per le metriche di cv
                          % Le possibilità sono:       'sig_mono', 'sig_singolo_diff', 'sig_doppio_diff'


IED = 5;                                            % IED espressa in mm
%% =========================================================================



%% ======================== Parametri filtraggio ===========================
tipo_filtro = "cheby2";
f_sample = 2048;                                    % Frequenza campionamento
f_taglio_basso = 20;                                % Frequenza minima del passabanda
f_taglio_alta = 400;                                % Frequenza massima del passabanda
f_notch = 50;                                       % Frequenza del notch
f_envelope = 4;                                     % Frequenza inviluppo
percH = 1.3;                                        % Percentuale frequenza alta
visualisation = "no";                               % Mostra grafici filtraggio

%% =========================================================================





%%  ========================= Import dati da file OTB e plot iniziali  =========================
% (codice docente con modifiche segnalate per evitare problemi con le celle)

FILTERSPEC = {'*.otb+','OTB+ files'; '*.otb','OTB file'; '*.zip', 'zip file'};
[FILENAME, PATHNAME] = uigetfile(FILTERSPEC,'titolo');

% Make new folder
mkdir('tmpopen');  
%cd('tempopen');

% Extract contents of tar file
untar([PATHNAME FILENAME],'tmpopen');
% Unzip([PATHNAME FILENAME]);
signals=dir(fullfile('tmpopen','*.sig')); %List folder contents and build full file name from parts
for nSig=1%:length(signals)
    PowerSupply{nSig}=3.3;
    abstracts{nSig}=[signals(nSig).name(1:end-4) '.xml'];
    abs = xml2struct(fullfile('.','tmpopen',abstracts{nSig}));

    %--- Modified section: Correcting the access to attributes ---
    if isfield(abs.Device, 'Attributes')
        Fsample{nSig}=str2num(abs.Device.Attributes.SampleFrequency);
        nChannel{nSig}=str2num(abs.Device.Attributes.DeviceTotalChannels);
        nADBit{nSig}=str2num(abs.Device.Attributes.ad_bits);
    end
    %--- End of modified section ---

    vett=zeros(1,nChannel{nSig});
    Gains{nSig}=vett;

    %--- Modified section: Ensure Adapter is handled as a cell array ---
    Adapter = abs.Device.Channels.Adapter;
    if ~iscell(Adapter)
        Adapter = {Adapter};
    end
    %--- End of modified section ---

    for nChild=1:length(Adapter)
        %--- Modified section: Correcting the access to gain and start index ---
        localGain{nSig}=str2num(Adapter{nChild}.Attributes.Gain);
        startIndex{nSig}=str2num(Adapter{nChild}.Attributes.ChannelStartIndex);
        %--- End of modified section ---

        Channel = Adapter{nChild}.Channel;
        if ~iscell(Channel)
            Channel = {Channel};
        end

        for nChan=1:length(Channel)
            ChannelAtt = Channel{nChan}.Attributes;
            idx=str2num(ChannelAtt.Index);
            Gains{nSig}(startIndex{nSig}+idx+1)=localGain{nSig};
        end
    end

    h=fopen(fullfile('tmpopen',signals(nSig).name),'r');
    data=fread(h,[nChannel{nSig} Inf],'short'); 
    fclose(h);

    processed(signals);

    Data{nSig}=data;
    figs{nSig}=figure;
    for nCh=1:nChannel{nSig}
       data(nCh,:)=data(nCh,:)*PowerSupply{nSig}/(2^nADBit{nSig})*1000/Gains{nSig}(nCh);
    end

    if mostra_plot_import
        MyPlotNormalized(figs{nSig},[1:length(data(1,:))]/Fsample{nSig},data);
        MyPlot(figure,[1:length(data(1,:))]/Fsample{nSig},data,0.5);
    end

end

rmdir('tmpopen','s');
%% =========================================================================



%% ========================= Filtraggio segnale =========================
data = data';
n_samples = size(data,1);
n_channel = size(data,2); 
sig_mono= zeros(length(data),n_channel);

% Filtraggio segnale
for i=1:n_channel
    sig_mono(:,i) = filter_general(data(:,i),tipo_filtro,f_sample,"fL",f_taglio_basso,"fH",f_taglio_alta,"fN",f_notch,"visualisation",visualisation);
end
%% =========================================================================



%% ========================= Creazione segnale singolo e doppio differenziale =========================

num_canali_singolo_diff = n_channel/2;
num_canali_doppio_diff = num_canali_singolo_diff/2;

% Inizializzazione della matrice dei segnali differenziali
sig_singolo_diff = zeros(n_samples, num_canali_singolo_diff);
sig_doppio_diff = zeros(n_samples, num_canali_doppio_diff);

% Popola matrice singolo differenziale
for i = 1:num_canali_singolo_diff
    canale_1 = 2 * i - 1;
    canale_2 = 2 * i;
    sig_singolo_diff(:, i) = sig_mono(:, canale_2) - sig_mono(:, canale_1);
end

% Popola matrice doppio differenziale
for i = 1:num_canali_doppio_diff
    canale_1 = 2 * i - 1;
    canale_2 = 2 * i;
    sig_doppio_diff(:, i) = sig_singolo_diff(:, canale_2) - sig_singolo_diff(:, canale_1);
end

% Plot di controllo mono, singolo e doppio
if mostra_plot_differenziali
    if is_bicipite
        subplot(3,1,1)
        plot(sig_mono(:,canale_plot))
        subtitle(['Segnale monopolare canale ' num2str(canale_plot)]);
    
        subplot(3,1,2)
        plot(sig_singolo_diff(:,canale_plot))
        subtitle(['Segnale singolo differenziale canale ' num2str(canale_plot)]);
    
        subplot(3,1,3)
        plot(sig_doppio_diff(:,canale_plot))
        subtitle(['Segnale doppio differenziale canale ' num2str(canale_plot)]);
    else
        subplot(3,1,1)
        plot(sig_mono(:,canale_plot+offset_canali))
        subtitle(['Segnale monopolare canale ' num2str(canale_plot)]);
    
        subplot(3,1,2)
        plot(sig_singolo_diff(:,canale_plot+offset_canali/2))       % Correzione valore avendo meno colonne
        subtitle(['Segnale singolo differenziale canale ' num2str(canale_plot)]);
    
        subplot(3,1,3)
        plot(sig_doppio_diff(:,canale_plot+offset_canali/4))        % Correzione valore avendo meno colonne
        subtitle(['Segnale doppio differenziale canale ' num2str(canale_plot)]);

    end
end
%% =========================================================================




%% ========================= Creazione fatigue plot=========================

% Segmentazione del segnale e calcolo valori per ogni segmento
passo = f_sample*risoluzione_calcolo;
numero_finestre = floor(n_samples/passo);

% Carica i segnali indicati all'inizio per effettuarne la successiva anlisi
eval(['segnale_metriche = ', nome_segnale_calcolo_metriche, ';']);
eval(['segnale_cv = ', nome_segnale_calcolo_cv, ';']);

if separa_muscoli
    rms = zeros(n_channel/2, numero_finestre);  % Diviso 2 perché intanto si lavora sempre e solo su metà dei sensori
    arv = zeros(n_channel/2,numero_finestre);
    mnf = zeros(n_channel/2,numero_finestre);
    mdf = zeros(n_channel/2,numero_finestre);
    cv = zeros(size(segnale_cv,2)/2,numero_finestre);
else
    rms = zeros(n_channel, numero_finestre);
    arv = zeros(n_channel,numero_finestre);
    mnf = zeros(n_channel,numero_finestre);
    mdf = zeros(n_channel,numero_finestre);
    cv = zeros(size(segnale_cv,2),numero_finestre);
end

for j = 1:numero_finestre
    start_idx = (j-1) * passo + 1;
    end_idx = start_idx + passo - 1;
    
    if end_idx > length(segnale_metriche)
        end_idx = length(segnale_metriche);
    end

    % Segmenta il segnale nella zona di interesse, utilizzando solo i
    % canali di interesse
    if separa_muscoli
        if is_bicipite
            segment_metrics = segnale_metriche(start_idx:end_idx, 1:n_channel/2);       
            segment_cv = segnale_cv(start_idx:end_idx, 1:size(segnale_cv,2)/2);
        else
            segment_metrics = segnale_metriche(start_idx:end_idx, n_channel/2+1:n_channel);
            segment_cv = segnale_cv(start_idx:end_idx, size(segnale_cv,2)/2+1:size(segnale_cv,2));
        end
    else
        segment_metrics = segnale_metriche(start_idx:end_idx, :);      
        segment_cv = segnale_cv(start_idx:end_idx, :);
    end

    % Calcola le varie metriche per la finedtra corretta
    [rms(:, j), arv(:, j), mnf(:, j), mdf(:, j), cv(:,j)] = FatiguePlot(segment_metrics, f_sample, IED, segment_cv);
end

%cv = sqrt(cv.*cv); % Ci sono dei valori negativi, DA CAPIRE PERCHÉ

asse_tempi = (0:numero_finestre-1) * risoluzione_calcolo;

% Fatigue plot normalizzati
if mostra_fatigue_normalizzati_singoli
    figure
    plot(asse_tempi,rms./rms(:,1))
    hold on
    plot(asse_tempi, arv./arv(:,1))
    hold on
    plot(asse_tempi, mnf./mnf(:,1))
    hold on 
    plot(asse_tempi, mdf./mdf(:,1))
    hold on
    plot(asse_tempi, cv/cv(1))
    legend('RMS', 'ARV', 'MNF', 'MDF', 'CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Calcoli andamenti medi tra i vari canali
mean_rms = mean(rms);   % Mean lavora sulla prima dimensione non unitaria
mean_arv = mean(arv);
mean_mnf = mean(mnf);
mean_mdf = mean(mdf);
mean_cv = mean(cv);           % Si lascia così per leggibilità del codice successivo, per come è calcolato cv è già un valore medio, avendo usato la stima multicanale

% Fatigue plot valori medi normalizzati
if mostra_fatigue_normalizzati__medi
    figure
    plot(asse_tempi,mean_rms/mean_rms(1))
    hold on
    plot(asse_tempi, mean_arv/mean_arv(1))
    hold on
    plot(asse_tempi, mean_mnf./mean_mnf(1))
    hold on 
    plot(asse_tempi, mean_mdf/mean_mdf(1))
    hold on
    plot(asse_tempi, mean_cv/mean_cv(1))
    
    legend('Mean RMS', 'Mean ARV', 'Mean MNF', 'Mean MDF', 'Mean CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Fatigue plot con subplot, valori non normalizzati
if mostra_fatigue_singoli
    figure 
    subplot(2,2,1)
    plot(asse_tempi,rms)
    subtitle("Valori RMS per canale")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')
    
    subplot(2,2,2)
    plot(asse_tempi,cv)
    subtitle("Valore medio CV")
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    
    subplot(2,2,3)
    plot(asse_tempi,mdf)
    subtitle("Valori MDF per canale")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    
    subplot(2,2,4)
    plot(asse_tempi,mnf)
    subtitle("Valori MNF per canale")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
end

if mostra_fatigue_medi   % Dovrebbero essere i grafici richiesti dal Mister
    figure 
    subplot(1,5,1)
    plot(asse_tempi,mean_rms)
    subtitle("Valore medio RMS")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')

    subplot(1,5,2)
    plot(asse_tempi,mean_arv)
    subtitle("Valore medio ARV")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')
    
    subplot(1,5,3)
    plot(asse_tempi,mean_cv)
    subtitle("Valore medio CV")
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    
    subplot(1,5,4)
    plot(asse_tempi,mean_mdf)
    subtitle("Valore medio MDF")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    
    subplot(1,5,5)
    plot(asse_tempi,mean_mnf)
    subtitle("Valore medio MNF")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
end
%% =========================================================================





%% ========================= Definizioni funzioni custom =========================

% Presa dal loro script
function []=MyPlotNormalized(fig,x,y)
    figure(fig);
    maximus=max(max(abs(y)));
    for ii=1:size(y,1)
        plot(x,y(ii,:)/2/maximus-ii);
        hold on
    end

end

% Presa dal loro script
function []=MyPlot(fig,x,y,shift)
    figure(fig);
    maximus=max(max(abs(y)));
    for ii=1:size(y,1)
        plot(x,y(ii,:)-ii*shift);
        hold on
    end
end