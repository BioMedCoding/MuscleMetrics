%% ========================== Initial notes ================================
% At the current state, the estimation of the CV appears to be incorrect; 
% there is one channel in particular that reports non-physiological values.



%% ========================== Initialization =============================
clear
close all
clc



%% ========================= General parameters ============================
%mostra_plot_import = true;
show_import_plot = true;
%mostra_plot_differenziali = true;
show_differential_plot = true;

%mostra_fatigue_normalizzati_singoli = true;
show_normalized_single_fatigue = true;
%mostra_fatigue_normalizzati__medi = true;
show_normalized_mean_fatigue = true;
%mostra_fatigue_singoli = true;
show_single_fatigue = true;
%mostra_fatigue_medi = true;
show_mean_fatigue = true;

%separa_muscoli = true;      % Se vero, va a passare alla funzione di analisi solamente la porzione di analisi definita dalla variabile is_bicipite
separate_muscles = true;      % If true, it will pass to the analysis function only the portion of analysis defined by the variable is_biceps.

is_biceps = true;         
% If true, only the first 8 channels of the matrix are passed to the analysis; otherwise, the last 8 channels are passed.
plot_channel = 1;            % Until resolved, it must be between 1 and 4 to avoid problems.
offset_channel = 8;          % Column separating the first and second signal (bi and tri)


%risoluzione_calcolo = 0.5;  % Risoluzione del fatigue plot, in secondi
fatigue_resolution = 0.5;  % Fatigue plot resolution in seconds

%nome_segnale_calcolo_metriche
metrics_sig_name = 'sig_singolo_diff';     % Nome del segnale che si vuole utilizzare nella parte dedicata ai fatigue plot per le metriche di ampiezza e frequenza
                                % Le possibilità sono: 'sig_mono', 'sig_singolo_diff', 'sig_doppio_diff'

%nome_segnale_calcolo_cv
cv_sig_name = 'sig_singolo_diff';   % Name of the signal to be used in the section dedicated to fatigue plots for CV metrics.
                          % Available possibilities are:       'sig_mono', 'sig_singolo_diff', 'sig_doppio_diff'


IED = 5;                                            % IED in mm
%% =========================================================================



%% ======================== Parametri filtraggio ===========================
filter_type = "cheby2";
f_sample = 2048;                                    % Sample frequency
f_cut_low = 20;                                     % Minimum frequency of the passband filter
f_cut_high = 400;                                   % Maximum frequency of the passband filter
f_notch = 50;                                       % Notch frequency
f_envelope = 4;                                     % Envelope frequency
percH = 1.3;                                        % High frequency percentage (to specify the filter caracteristic)
visualisation = "no";                               % Show filtering plots

%% =========================================================================





%%  ========================= Import data from OTB file and initial plots  =========================
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

    if show_import_plot
        MyPlotNormalized(figs{nSig},[1:length(data(1,:))]/Fsample{nSig},data);
        MyPlot(figure,[1:length(data(1,:))]/Fsample{nSig},data,0.5);
    end

end

rmdir('tmpopen','s');
%% =========================================================================



%% ========================= Signal filtering =========================
data = data';
n_samples = size(data,1);
n_channel = size(data,2); 
sig_mono= zeros(length(data),n_channel);

% Signal filtering
for i=1:n_channel
    sig_mono(:,i) = filter_general(data(:,i),filter_type,f_sample,"fL",f_cut_low,"fH",f_cut_high,"fN",f_notch,"visualisation",visualisation);
end
%% =========================================================================



%% ========================= Creation of single and double differential signals =========================

n_channel_single_diff = n_channel/2;
n_channel_double_diff = n_channel_single_diff/2;

% Initialization of the differential signals matrix
sig_single_diff = zeros(n_samples, n_channel_single_diff);
sig_double_diff = zeros(n_samples, n_channel_double_diff);

% Populate single differential matrix
for i = 1:n_channel_single_diff
    channel_1 = 2 * i - 1;
    channel_2 = 2 * i;
    sig_single_diff(:, i) = sig_mono(:, channel_2) - sig_mono(:, channel_1);
end

% Populate double differential matrix
for i = 1:n_channel_double_diff
    channel_1 = 2 * i - 1;
    channel_2 = 2 * i;
    sig_double_diff(:, i) = sig_single_diff(:, channel_2) - sig_single_diff(:, channel_1);
end

% Control plot for mono, single, and double differential signals
if show_differential_plot
    if is_biceps
        subplot(3,1,1)
        plot(sig_mono(:,plot_channel))
        subtitle(['Segnale monopolare canale ' num2str(plot_channel)]);
    
        subplot(3,1,2)
        plot(sig_single_diff(:,plot_channel))
        subtitle(['Segnale singolo differenziale canale ' num2str(plot_channel)]);
    
        subplot(3,1,3)
        plot(sig_double_diff(:,plot_channel))
        subtitle(['Segnale doppio differenziale canale ' num2str(plot_channel)]);
    else
        subplot(3,1,1)
        plot(sig_mono(:,plot_channel+offset_channel))
        subtitle(['Segnale monopolare canale ' num2str(plot_channel)]);
    
        subplot(3,1,2)
        plot(sig_single_diff(:,plot_channel+offset_channel/2))       % Correzione valore avendo meno colonne
        subtitle(['Segnale singolo differenziale canale ' num2str(plot_channel)]);
    
        subplot(3,1,3)
        plot(sig_double_diff(:,plot_channel+offset_channel/4))        % Correzione valore avendo meno colonne
        subtitle(['Segnale doppio differenziale canale ' num2str(plot_channel)]);

    end
end
%% =========================================================================




%% ========================= Fatigue plot creation =========================

% Signal segmentation variables
step = f_sample*fatigue_resolution;
n_window = floor(n_samples/step);

% Load the signals indicated at the beginning for subsequent analysis
eval(['segnale_metriche = ', metrics_sig_name, ';']);
eval(['segnale_cv = ', cv_sig_name, ';']);

if separate_muscles
    rms = zeros(n_channel/2, n_window);  % Divided by 2 because we are always working with only half of the sensors
    arv = zeros(n_channel/2,n_window);
    mnf = zeros(n_channel/2,n_window);
    mdf = zeros(n_channel/2,n_window);
    cv = zeros(size(segnale_cv,2)/2,n_window);
else
    rms = zeros(n_channel, n_window);
    arv = zeros(n_channel,n_window);
    mnf = zeros(n_channel,n_window);
    mdf = zeros(n_channel,n_window);
    cv = zeros(size(segnale_cv,2),n_window);
end

for j = 1:n_window
    start_idx = (j-1) * step + 1;
    end_idx = start_idx + step - 1;
    
    if end_idx > length(segnale_metriche)
        end_idx = length(segnale_metriche);
    end

    % Segment the signal in the region of interest, using only the channels of interest.
    if separate_muscles
        if is_biceps
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

    % Calculate the various metrics for the correct window.
    [rms(:, j), arv(:, j), mnf(:, j), mdf(:, j), cv(:,j)] = FatiguePlot(segment_metrics, f_sample, IED, segment_cv);
end

%cv = sqrt(cv.*cv); % Ci sono dei valori negativi, DA CAPIRE PERCHÉ

time_axis = (0:n_window-1) * fatigue_resolution;

% Normalized fatigue plot
if show_normalized_single_fatigue
    figure
    plot(time_axis,rms./rms(:,1))
    hold on
    plot(time_axis, arv./arv(:,1))
    hold on
    plot(time_axis, mnf./mnf(:,1))
    hold on 
    plot(time_axis, mdf./mdf(:,1))
    hold on
    plot(time_axis, cv/cv(1))
    legend('RMS', 'ARV', 'MNF', 'MDF', 'CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Calculate mean values across the various channels
mean_rms = mean(rms);   % Mean works on the first non unitary dimension
mean_arv = mean(arv);
mean_mnf = mean(mnf);
mean_mdf = mean(mdf);
mean_cv = mean(cv);           

% Mean normalized fatigue plot
if show_normalized_mean_fatigue
    figure
    plot(time_axis,mean_rms/mean_rms(1))
    hold on
    plot(time_axis, mean_arv/mean_arv(1))
    hold on
    plot(time_axis, mean_mnf./mean_mnf(1))
    hold on 
    plot(time_axis, mean_mdf/mean_mdf(1))
    hold on
    plot(time_axis, mean_cv/mean_cv(1))
    
    legend('Mean RMS', 'Mean ARV', 'Mean MNF', 'Mean MDF', 'Mean CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Fatigue plot with subplot, non-normalized values
if show_single_fatigue
    figure 
    subplot(2,2,1)
    plot(time_axis,rms)
    subtitle("Valori RMS per canale")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')
    
    subplot(2,2,2)
    plot(time_axis,cv)
    subtitle("Valore medio CV")
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    
    subplot(2,2,3)
    plot(time_axis,mdf)
    subtitle("Valori MDF per canale")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    
    subplot(2,2,4)
    plot(time_axis,mnf)
    subtitle("Valori MNF per canale")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
end

if show_mean_fatigue   % Dovrebbero essere i grafici richiesti dal Mister
    figure 
    subplot(1,5,1)
    plot(time_axis,mean_rms)
    subtitle("Valore medio RMS")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')

    subplot(1,5,2)
    plot(time_axis,mean_arv)
    subtitle("Valore medio ARV")
    xlabel('Time [s]')
    ylabel('Amplitude [mv?]')
    
    subplot(1,5,3)
    plot(time_axis,mean_cv)
    subtitle("Valore medio CV")
    xlabel('Time [s]')
    ylabel('Speed [m/s]')
    
    subplot(1,5,4)
    plot(time_axis,mean_mdf)
    subtitle("Valore medio MDF")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    
    subplot(1,5,5)
    plot(time_axis,mean_mnf)
    subtitle("Valore medio MNF")
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
end
%% =========================================================================





%% ========================= Custom function definition =========================

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