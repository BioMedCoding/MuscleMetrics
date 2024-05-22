function aa=daf(x,ntlag,fc,T)
% DAF - Discrete Ambiguity Function (\tau,\theta), both ordered with
% positive value first and then the negative ones flipped
% To order the result: aa=fftshift(aa);
% Usage:   aa=daf(x,ntlag,fc,T);
% 
% Inputs:  x      input signal
%          ntlag  number of time lag for the autocorrelation sequence estimate
%          fc     sampling rate in Hz
%          T      duration of the signal in seconds
% Output:  aa     matrix containing the DAF


y=act(x,ntlag); % (\tau,t)
% direct Fourier transform t -> \theta
aa=(1/fc)*(1/T)*fft(y,[],2); % (\tau,\theta)
