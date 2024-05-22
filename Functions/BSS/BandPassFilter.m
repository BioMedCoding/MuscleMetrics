function sigfil=BandPassFilter(XX,fL,fH,fc)

% Function performing a band pass filter of the signals contained in the
% matrix XX
%
% Input:
% XX: matrix containing on each row signal
% fL, fH: cutoff frequencies (f1: low frequency, f2: high frequency)
% fc: sampling frequency
%
% Output:
% sigfil: filtered signal
f2=fL;
f1=.75*fL;
Wp=2*f2/fc;
Ws=2*f1/fc;
Rp=1;
Rs=20;
[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
[BH,AH] = cheby2(n,Rs,Ws,'high');

f2=fH;
f1=1.1*fH;
Wp=2*f2/fc;
Ws=2*f1/fc;
Rp=1;
Rs=20;[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
[BL,AL] = cheby2(n,Rs,Ws);

for ch = 1:size(XX,1)
    sigfil(ch,:) = filtfilt(BH,AH,XX(ch,:));
    sigfil(ch,:) = filtfilt(BL,AL,sigfil(ch,:));
end