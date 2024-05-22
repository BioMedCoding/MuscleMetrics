% Applicazione dell'ICA alla rimozione di artefatti.
% Vengono simulate 3 sorgenti: 2 rumori (segnali di nostro interesse) e una
% componente di bassa frequenza (artefatto da togliere).
% Si usa l'ICA per separare le 3 componenti, identificare l'artefatto e
% rimuoverlo.

close all;
clear
clc

%% mixing matrix
M=[-1 2 3; 2 -2 1; 2 -1 -2];

%% sorgenti
N=1000;
fs=100;
t=0:1/fs:N/fs;
t=t-mean(t);

%% generare 2 sorgenti di N campioni con distribuzione uniforme
x1 = randn(1, N);
x3 = randn(1,N);

%% Sorgente di bassa frequenza (artefatto da rimuovere)
sigma = 1;
x2 = exp(-t.^2/2/sigma^2);
% Vettore delle sorgenti (mettere in un vettore le 3 sorgenti simulate)
x = vertcat(x1, x2, x3);

% normalizzare le sorgenti usando la funzione zscore
for i=1:3
    x(i,:) = zscore(x(i,:));
end

% mixtures
y=M*x;

% ICA
[icasig, A, W] = fastica(y);

% Calcolare le IC usando le miscele e la matrice di de-mixing
ic = W * y;

% controllare che le IC calcolate nella riga precedente siano uguali alle
% icasig
for i=1:3
    subplot(3,1,i);
    plot(icasig(i,:));
    hold on;plot(ic(i,:));
    title(['IC' num2str(i)]);
end

% Ricostruzione delle miscele usando le IC e la matrice di mixing
yest= A * icasig;

% test che le miscele ricostruite siano davvero uguali alle miscele
% iniziali
figure;
for i=1:3
    subplot(3,1,i);
    plot(y(i,:));
    hold on
    plot(yest(i,:));
    title(['mixture ' num2str(i)]);
end

% rimozione della componente di bassa frequenza
% Trovare l'indice index_IC della componente di bassa frequenza 
%index_IC=...;
[~, index_IC] = max(abs(corr(x(2, :)', icasig')));

% vettore delle sorgenti buone (sottraggo index_IC dal vettore 1:3)
vet_good_IC=setdiff(1:3,index_IC);

% ricostruzione delle miscele usando solo le sorgenti buone
%y_artifact_removed=...;
icasig_no_artifact = icasig(vet_good_IC, :);
A_no_artifact = A(:, vet_good_IC);
y_artifact_removed = A_no_artifact * icasig_no_artifact;

% Dati senza l'artefatto di bassa frequenza
figure;
for i=1:3
    subplot(3,1,i);
    plot(y_artifact_removed(i,:));
    title(['mixture ' num2str(i) ' without the artifact']);
end

