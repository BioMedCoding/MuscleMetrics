function w=dcw(x,fc,ntlag,sig,T)
% DCW - Discrete Choi-Williams time-frequency transform
% Usage:   w=dcw(x,fc,ntlag,sig,T);
% 
% Inputs:  x      input signal
%          ntlag  number of time lag for the autocorrelation sequence estimate
%          fc     sampling rate in Hz
%          sig    value of the selectivity parameter (>0 and <100)
%          T      time support of the input signal in seconds
% Output:  w      matrix containing the Choi-Williams transform

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised

aa=1/((1/fc)*(1/T))*daf(x,ntlag,fc,T);
[rowaa, colaa]=size(aa);
kcw=fftshift(kercw(sig,rowaa,colaa,fc));
ww=kcw.*aa;clear aa kcw
% inverse transform \theta -> t
w=ifft(ww,[],2);clear ww
% direct transform \tau -> f
www=1/(ntlag*2/fc)*(2/fc)*fft(w);clear w
w=www;clear www
[roww , ~]=size(w);
w=w(1:roww/2+1,:);
