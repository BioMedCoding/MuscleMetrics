function [rms, arv, mnf, mdf, cv_total, cv_overall] = FatiguePlot(sig_valori, fs, IED, sig_cv)

% input
% sig: matrix N x M, with N samples and N channels


% Preinizializza variabili
num_channels = size(sig_valori, 2); % Numero di canali, se rispettata sintassi
rms = zeros(1, num_channels);
arv = zeros(1, num_channels);
mnf = zeros(1, num_channels);
mdf = zeros(1, num_channels);
cv = zeros(1, size(sig_cv,2)-1); % -1 aggiunto per non avere l'ultimo valore di 0

% Iera calcolo sui singoli canali
for ch = 1:num_channels

    current_sig = sig_valori(:, ch);
    
    % Metriche di ampiezza
    rms(ch) = sqrt(mean(current_sig.^2));   % RMS
    arv(ch) = sum(abs(current_sig)) / length(current_sig);  % ARV
    
    % Metriche di frequenza
    [Pxx, f] = pwelch(current_sig - mean(current_sig), tukeywin(length(current_sig), .1), 0, [], fs);
    mnf(ch) = sum(f .* Pxx) / sum(Pxx); % Mean Frequency (MNF)
    
    cumulative_psd = cumsum(Pxx);
    median_value = cumulative_psd(end) / 2;
    median_index = find(cumulative_psd >= median_value, 1);
    mdf(ch) = f(median_index); % Median Frequency (MDF)
    
   

end

clear ch num_channels

 % Velocità di conduzione
num_channels = size(sig_cv, 2);

for ch = 1:(num_channels-1)
        s1 = sig_cv(:,ch);
        s2 = sig_cv(:,ch+1);
        [xc, del] = xcorr(s2, s1, 6);
        [mm, I] = max(xc);
        start = del(I);
        d = delay(real(fft(s2)), imag(fft(s2)), real(fft(s1)), imag(fft(s1)), start);
        cv(ch) = abs(IED / (d * 1000) * fs);
end

cv_total = cv;

%Per adattarlo alla funzione
IED = IED/1000;
cv_overall = mle_CV_est(sig_cv(:,:)', IED, fs);

end
