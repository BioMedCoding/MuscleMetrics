%% Inizializzazione
clear
close all
clc

%% ======================== Parametri filtraggio ===========================
tipo_filtro = "cheby2";
f_sample = 2048;                                    % Frequenza campionamento
f_taglio_basso = 20;                                % Frequenza minima del passabanda
f_taglio_alta = 400;                                % Frequenza massima del passabanda
f_notch = 50;                                       % Frequenza del notch
f_envelope = 4;                                     % Frequenza inviluppo
percH = 1.3;                                        % Percentuale frequenza alta
visualisation = "no";                               % Mostra grafici filtraggio

IED = 5;                                            % IED espressa in mm
%% =========================================================================


%% Import dati da file OTB e plot iniziali

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

    MyPlotNormalized(figs{nSig},[1:length(data(1,:))]/Fsample{nSig},data);
    MyPlot(figure,[1:length(data(1,:))]/Fsample{nSig},data,0.5);

end

rmdir('tmpopen','s');

% % theFiles = dir;
% % for k = 1 : length(theFiles)
% %   baseFileName = theFiles(k).name;
% %   fprintf(1, 'Now deleting %s\n', baseFileName);
% %   delete(baseFileName);
% % end
% % 
% % cd ..;
% % rmdir('tempopen','s');

%% Filtraggio segnale 
data = data';
n_channel = length(data(1,:));
sig_filt= zeros(length(data),n_channel);

% Filtraggio segnale
for i=1:n_channel
    sig_filt(:,i) = filter_general(data(:,i),tipo_filtro,f_sample,"fL",f_taglio_basso,"fH",f_taglio_alta,"fN",f_notch,"visualisation",visualisation);
end

%% Creazione fatigue plot
s1 = -2;
s2 = -2;
[rms,arv,mnf,mdf,cv] = FatiguePlot(data,s1,s2,f_sample,IED);

%% Definizioni funzioni custom

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