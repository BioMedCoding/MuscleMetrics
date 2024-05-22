clear;
close all
clc


fc = 1000; % Frequenza di campionamento supposta
nsamples=100;
fL=[2 100];
fH=[300 350];


%% generare 2 segnali casuali, filtrando un rumore Gaussiano nelle bande 2-300 Hz e 100-350 Hz
    %input_signal(1,:) = ...;
    %input_signal(2,:) = ...;

    % Seguendo il genio
    noise1 = randn(1,nsamples);
    noise2 = randn(1, nsamples);

    input_signal(1,:) = BandPassFilter(noise1, fL(1), fH(1), fc);
    %input_signal(1,:) = BandPassFilter(input_signal(1,:), fL(2), fH(2), fc);

    input_signal(2,:) = BandPassFilter(noise2, fL(2), fH(2), fc);
    %input_signal(2,:) = BandPassFilter(input_signal(2,:), fL(2), fH(2), fc);

%% faccio in modo che i 2 segnali siano ortogonali    
    % Basato su ortonormalizzazione di Schmitt
input_signal(1,:)=input_signal(1,:)-(input_signal(1,:)*input_signal(2,:)')/length(input_signal)*input_signal(2,:);

%% definisco il modello di mixing e il rumore
N=2;
noise=0.01*randn(N,nsamples);% rumore

%% Creazione matrice di mixing 
% M=randn(N,2); % matrice di mixing casuale
M=[.25 .25;-0.1 1]; % matrice di mixing scelta dall'utente

%% Miscele
noisy_input=M*input_signal+noise; 

[coeff,s,l]=pca(noisy_input'); % Usare la funzione Matlab pca per separare le sorgenti
    % Calcoli matrice di autocorrelazione disponendo correttamente le
% matrici, c'è già la funzione di Matlab
output = s';           % Calcolare le componenti principali in funzione del tempo
output_manuale = coeff*noisy_input;

figure
plot(output(1,:))
hold on 
plot(output_manuale(1,:))

%% normalizzazione delle componenti principali
for i=1:2
    output(i,:)=output(i,:)/std(output(i,:));
end

%% Figure
figure;
subplot(3,2,1);
plot(input_signal(1,:));aa=axis;%
title('Original Component 1')
subplot(3,2,2);
plot(input_signal(2,:), 'r');bb=axis;    % plot B
title('Original Component 2')

subplot(3,2,3);
for i=1:N
    plot(noisy_input(i,:)/range(noisy_input(:))+i);
    hold on
end      % plot mixing 1
title('Input Data Set')

subplot(3,2,4)
plot(noisy_input(1,:),noisy_input(2,:),'k.');
hold on
plot([0 coeff(1,1)],[0 coeff(1,2)],'r')
plot([0 coeff(2,1)],[0 coeff(2,2)],'r')
axis equal
title('Scatter plot and directions of PCs')

%% Mettere in ordine le componenti in base alla correlazione con le sorgenti
% input_signal
RR = corr(input_signal', output'); % coefficiente di correlazione - completato dal genio
I=1:2;
[~,J(1)]=max(abs(RR(1,:)));
[~,J(2)]=max(abs(RR(2,:)));

for i=1:2
    Y(i,:)=round(RR(I(i),J(i)))*output(J(i),:);
end
subplot(3,2,5);
plot(Y(1,:),'b'); axis(aa)
title('Retreived Component-1') %Component 1

subplot(3,2,6);
plot(Y(2,:),'r');axis(bb)
title('Retreived Component-2') %Component 2

linkaxes([subplot(3,2,1), subplot(3,2,2), subplot(3,2,5), subplot(3,2,6)])