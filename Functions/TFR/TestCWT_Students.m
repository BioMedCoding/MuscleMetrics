% Test della routine per la trasformata wavelet continua
clear;
close all
clc

%% parametro che rientra nella definizione della mother wavelet (standard deviation)
global sd

%% Percentuale di energia che si accetta di scartare
perc_discard=5;
fsamp=2048;
sd=28/fsamp/10;
Nlevels=500;

% Carico un segnale EMG SD simulato e ne prendo una porzione
load('EMG_4WT');
sig=SD(2,2*fsamp+1:3*fsamp);

% taglio l'inizio e la fine del segnale per togliere gli effetti di bordo
sig=sig.*tukeywin(length(sig),.1)';

% definizione del vettore tempo
N=length(sig);
t=[1:N]/fsamp;

% la media viene tolta (comunque non potrei ricostruirla con la trasformata wavelet)
sig=sig-mean(sig);

% calcolo della PSD (periodogramma semplice), per definire come filtrare il
% segnale in modo che le ondine considerate possano approssimarlo bene
[pxx,fr] = pwelch(...);
    
% Definizione della Mother Wavelet 

f=@(t) ...;% Derivata di una Gaussiana con standard deviation sd
    
figure;
vt=length(t)/2+1-20:length(t)/2+20;plot(t(vt),f(t(vt)));
title('Mother wavelet');
axis off

% PSD della mother wavelet
[pxx_f,f_f] = pwelch(f(t),...);
    % Frequenza media della mother wavelet (serve per definire la relazione tra scala e frequenza)
f_0=mean(f_f.*pxx_f)/mean(pxx_f);

% Filtro passa-banda in modo da buttare via una percentuale di
% energia pari a perc_discard (definito in precedenza). Vengono eliminate
% così delle basse e alte frequenze che non potrei rappresentare con le
% ondine considerate

% 1. Usare la funzione cumsum per fare l'integrale della PSD
vv=cumsum(pxx)/sum(pxx);

% 2. Definire la frequenza minima imponendo di escludere una percentuale di
% potenza pari a perc_discard/2
fmin=...

% 3. Definire la frequenza massima imponendo di escludere una percentuale di
% potenza pari a perc_discard/2
fmax=...;
    % figura che mostra la PSD del segnale e dell'ondina. Variare sd finché la

% PSD dell'ondina "assomiglia" a quella del segnale
figure;subplot(2,1,2);plot(fr, pxx/max(pxx));hold on;plot(f_f,pxx_f/max(pxx_f));xlabel('Frequency (Hz)');title('PSDs of the signal and of the mother wavelet (and minimum and maximum frequency considered)');plot([fmin fmin],[0 1],'k--');plot([fmax fmax],[0 1],'k--')
subplot(2,1,1);plot(t,sig)

% filtro passabanda che tiene le frequenze comprese tra fmin e fmax
% 1. filtro passa-alto (usare ad esempio la funzione cheby2 per definire i pesi del filtro)
...
    % Filtro a doppia passata per eliminare la distorsione di fase
sig = filtfilt(b,a,sig')';
% 2. filtro passa-basso
...
    sig = filtfilt(b,a,sig')';
% Mostriamo il segnale filtrato
hold on;plot(t,sig);
legend('original','filtered');
xlabel('Time (s)');
title('Signal to be processed')

% Definizione delle scale minima e massima (in funzione della frequenza
% media della mother wavelet e del range di frequenza del segnale) 
scale_min=f_0/fmax;scale_max=f_0/fmin;

% Scegliere una distribuzione (ad esempio lineare) degli
% scalamenti delle ondine
yscale=linspace(...);
    
% Trasformata wavelet
[scalo, cwt] = CWTW_LM(sig,yscale,f,t);

% Grafico del risultato
figure;subplot(2,1,1);
plot(t,sig);title('Signal');
xlabel('Time (s)');
subplot(2,1,2);
image((flipud(scalo')-min(scalo(:)))/range(scalo(:))*256)
colormap(gray(255));
title('CWT');
xlabel('Time samples');
ylabel('Scales');

% ricostruzione del segnale tramite la trasformata inversa
y = ICWTW_LM(cwt, yscale,f,t);
figure;

% Mostro il segnale ricostruito e indico l'errore
plot(y,'b');
hold on;
plot(sig,'r');
title(['Signal and reconstruction - error = ' num2str(round(std(x(:)-sig(:))/std(sig)*10000)/100) '%']);

