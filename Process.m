%% ========================== Initial notes ================================
% At the current state, the estimation of the CV appears to be incorrect; 
% there is one channel in particular that reports non-physiological values.

% Trasposizione matrice segnale
% Dio cane la CV continua a non funzionare né con calcolo diretto né con
% calcolo complessivo


%% ========================== Initialization =============================
clear
close all
clc



%% ========================= General parameters ============================

show_import_plot = true;
show_differential_plot = true;

show_normalized_single_fatigue = true;
show_normalized_mean_fatigue = true;
show_single_fatigue = true;
show_mean_fatigue = true;

separate_muscles = true;              % If true, it will pass to the analysis function only the portion of analysis defined by the variable is_biceps.

is_biceps = false;                     % If true, only the first 8 channels of the matrix are passed to the analysis; otherwise, the last 8 channels are passed.
plot_channel = 1;                     % Until resolved, it must be between 1 and 4 to avoid problems.
muscle_channel = 8;                   % Column separating the first and second signal (bi and tri)

fatigue_resolution = 0.5;             % Fatigue plot resolution in seconds

metrics_sig_name = 'single_diff';     % Nome del segnale che si vuole utilizzare nella parte dedicata ai fatigue plot per le metriche di ampiezza e frequenza
                                      % Available possibilities are: 'mono', 'single_diff', 'double_diff'

cv_sig_name = 'single_diff';          % Name of the signal to be used in the section dedicated to fatigue plots for CV metrics.
                                      % Available possibilities are:       'mono', 'single_diff', 'double_diff'

IED = 5;                                            % IED in mm

remove_outliers = true;
%% =========================================================================



%% ======================== Filter parameters ===========================
filter_type = "cheby2";
f_sample = 2048;                       % Sample frequency
f_cut_low = 20;                        % Minimum frequency of the passband filter
f_cut_high = 400;                      % Maximum frequency of the passband filter
f_notch = 50;                          % Notch frequency
percH = 1.3;                           % High frequency percentage (to specify the filter caracteristic)
visualisation = "no";                  % Show filtering plots

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
n_channel_total = size(data,2);

bic = data(:, 1:muscle_channel);
tri = data(:, muscle_channel+1:n_channel_total);

bic_mono = zeros(n_samples, muscle_channel);
tri_mono = zeros(n_samples, muscle_channel);

% Signal filtering
for i=1:muscle_channel
    bic_mono(:,i) = filter_general(bic(:,i),filter_type,f_sample,"fL",f_cut_low,"fH",f_cut_high,"fN",f_notch,"visualisation",visualisation);
    tri_mono(:,i) = filter_general(tri(:,i),filter_type,f_sample,"fL",f_cut_low,"fH",f_cut_high,"fN",f_notch,"visualisation",visualisation);
end

%% =========================================================================



%% ========================= Creation of single and double differential signals =========================

n_channel_single_diff = size(bic_mono,2)-1;      % I due set intanto sono uguali
n_channel_double_diff = n_channel_single_diff-1;

% Create the single differential matrix
bic_single_diff = diff(bic_mono')';  % Doppia trasposizione per adattare a funzione e tornare ad originale
tri_single_diff = diff(tri_mono')';

% Create the double differential matrix
bic_double_diff = diff(bic_single_diff')';
tri_double_diff = diff(tri_single_diff')';

% Control plot for mono, single, and double differential signals
if show_differential_plot
    if is_biceps
        subplot(3,1,1)
        plot(bic_mono(:,plot_channel))
        subtitle(['Monopolar signal - channel ' num2str(plot_channel)]);
    
        subplot(3,1,2)
        plot(bic_single_diff(:,plot_channel))
        subtitle(['Single differential signal - channel ' num2str(plot_channel)]);
    
        subplot(3,1,3)
        plot(bic_double_diff(:,plot_channel))
        subtitle(['Double differential signal - channel ' num2str(plot_channel)]);

        linkaxes([subplot(3,1,1), subplot(3,1,2), subplot(3,1,3)], 'x')

    else
        subplot(3,1,1)
        plot(tri_mono(:,plot_channel))
        subtitle(['Monopolar signal - channel ' num2str(plot_channel)]);
    
        subplot(3,1,2)
        plot(tri_single_diff(:,plot_channel))       % Correzione valore avendo meno colonne
        subtitle(['Single differential signal - channel ' num2str(plot_channel)]);
    
        subplot(3,1,3)
        plot(tri_double_diff(:,plot_channel))        % Correzione valore avendo meno colonne
        subtitle(['Double differential signal - channel ' num2str(plot_channel)]);

        linkaxes([subplot(3,1,1), subplot(3,1,2), subplot(3,1,3)], 'x')


    end
end
%% =========================================================================




%% ========================= Fatigue plot creation =========================

% Signal segmentation variables
step = f_sample*fatigue_resolution;
n_window = floor(n_samples/step);

% Load the signals indicated at the beginning for subsequent analysis
if is_biceps
    eval(['metrics_signal = ','bic_', metrics_sig_name, ';']);
    eval(['cv_signal = ','bic_', cv_sig_name, ';']);
else
    eval(['metrics_signal = ','tri_', metrics_sig_name, ';']);
    eval(['cv_signal = ','tri_', cv_sig_name, ';']);
end

n_channel_metrics = size(metrics_signal,2);
n_channel_cv = size(cv_signal,2);

rms = zeros(n_channel_metrics, n_window);
arv = zeros(n_channel_metrics,n_window);
mnf = zeros(n_channel_metrics,n_window);
mdf = zeros(n_channel_metrics,n_window);
cv = zeros(n_channel_cv-1,n_window);


for j = 1:n_window
    start_idx = (j-1) * step + 1;
    end_idx = start_idx + step - 1;
    
    if end_idx > length(metrics_signal)
        end_idx = length(metrics_signal);
    end

    % Segment the signal in the region of interest, using only the channels of interest.
    segment_metrics = metrics_signal(start_idx:end_idx, :);       
    segment_cv = cv_signal(start_idx:end_idx, :);
   

    % Calculate the various metrics for the correct window.
    [rms(:, j), arv(:, j), mnf(:, j), mdf(:, j), cv(:,j)] = FatiguePlot(segment_metrics, f_sample, IED, segment_cv);
end

%cv = sqrt(cv.*cv); % Ci sono dei valori negativi, DA CAPIRE PERCHÉ
if remove_outliers
    cv_cleaned = excludeOutliers(cv,1); 
    mean_cv = mean(cv);
    mean_cv_cleaned = mean(cv_cleaned);

    figure
    hold on
    title('Comparison in mean cv values original vs outlier_cleaned')
    plot(mean_cv)
    plot(mean_cv_cleaned)
    legend('Original', 'Outliers_cleaned')

    cv = cv_cleaned;
end



time_axis = (0:n_window-1) * fatigue_resolution;

% Normalized fatigue plot
if show_normalized_single_fatigue
    
    figure
    hold on
    title('Normalized fatigue plot, single channels')

    plot(time_axis,rms./rms(:,1))
    plot(time_axis, arv./arv(:,1))
    plot(time_axis, mnf./mnf(:,1))
    plot(time_axis, mdf./mdf(:,1))
    plot(time_axis, cv/cv(1))
    legend('RMS', 'ARV', 'MNF', 'MDF', 'CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Exclude outliers values from the data

% Calculate mean values across the various channels
mean_rms = mean(rms);   % Mean works on the first non unitary dimension
mean_arv = mean(arv);
mean_mnf = mean(mnf);
mean_mdf = mean(mdf);
mean_cv = mean(cv); 
  

% Mean normalized fatigue plot
if show_normalized_mean_fatigue
    
    figure
    hold on
    title('Normalized fatigue plot, mean values')

    plot(time_axis,mean_rms/mean_rms(1))
    plot(time_axis, mean_arv/mean_arv(1))
    plot(time_axis, mean_mnf./mean_mnf(1))
    plot(time_axis, mean_mdf/mean_mdf(1))
    plot(time_axis, mean_cv/mean_cv(1))
    
    legend('Mean RMS', 'Mean ARV', 'Mean MNF', 'Mean MDF', 'Mean CV')
    xlabel('Time [s]')
    ylabel('Normalized unit')
end

% Fatigue plot with subplot, non-normalized values
if show_single_fatigue
    
    figure
    sgtitle('Fatigue plot, single channels')
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
    sgtitle('Fatigue plot, mean values')
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

function clean_matrix = excludeOutliers(original_matrix, dimension)
    
        % dimension 1 = columns, dimensione 2 = rows


        % Calcolo valore mediano 
        median_values = median(original_matrix, dimension);

        % Calcolo deviazione assoluta dalla mediana
        threshold = mad(original_matrix,dimension);

        lower_limit = max(median_values - 3*threshold,0);
        higher_limit = min(median_values +3*threshold,10);

        if dimension == 1
            direction_elimination = 2;
        else
            direction_elimination = 1;
        end

        % Identificare le righe che non contengono outlier
        rowsToKeep = all(original_matrix >= lower_limit & original_matrix <= higher_limit,direction_elimination);
    
        % Creare la matrice pulita escludendo le righe con outlier
        clean_matrix = original_matrix(rowsToKeep, :);

 end
    