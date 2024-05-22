function w=wig(x,fc,ntlag)
% WIG - Wigner-Ville time-frequency transform
% Usage:   w=wig(x,fc,ntlag);
% 
% Inputs:  x      input signal
%          ntlag  number of time lag for the autocorrelation sequence estimate
%          fc     sampling rate in Hz
% Output:  w      matrix containing the Wigner-Ville transform

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised

y=act(x,ntlag);
[rowy,~]=size(y);
w=1/(ntlag*2/fc)*(2/fc)*fft(y);
w=w(1:rowy/2+1,:);

