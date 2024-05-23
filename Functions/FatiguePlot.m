function [rms, arv, mnf, mdf, cv] = FatiguePlot(sig, s1, s2, fs, IED)

% input
% sig: matrix N x M, with N samples and N channels

% Preinizializza variabili
num_channels = size(sig, 2); % Numero di canali, se rispettata sintassi
rms = zeros(1, num_channels);
arv = zeros(1, num_channels);
mnf = zeros(1, num_channels);
mdf = zeros(1, num_channels);
cv = zeros(1, num_channels);

% Iera calcolo sui singoli canali
for ch = 1:num_channels

    current_sig = sig(:, ch);
    
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
    
    % Velocità di conduzione
    [xc, del] = xcorr(s2, s1, 10);
    [mm, I] = max(xc);
    start = del(I);
    d = delay(real(fft(s2)), imag(fft(s2)), real(fft(s1)), imag(fft(s1)), start);
    cv(ch) = IED / (d * 1000) * fs;
end

end
