function y=act(x,ntlag)
%act Instantaneous Autocorrelation Function 
% (\tau,t) with the positive delays \tau at the beginning and then the
% negative ones flipped (ready for fft)
% To get the solution with \tau ordered (first negative, then positive
% values) apply the following: y=fftshift(y,1);
%	y=act(x,ntlag) estimates the instantaneous autocorrelation function
%	of an input signal x. The number of time lags considered to estimate
%	the instantaneous autocorrelation function is equal to ntlag. It
%	corresponds to the number of rows of the output matrix unless
%	ntlag is lower than the number of samples of the input signal.
%	Each row of y contains the products x[n+k/2]*conj(x[n-k/2]) for a
%	value of k (shifts of half a sample are obtained shifting the signal in
%	the frequency domain), where n is related to the time axis and k is related
%	to the time lag axis. If the number of points of the
%	input vector is smaller than ntlag, the number of time lags
%	is limited to the number of samples of the input signal.

%	input variable(s):	x	input signal 
%				ntlag	number of time lags

%	output variable(s):	y	instantaneous autocorrelation function

%	author(s):	G. Gagliati
%			P. Bonato, 4-22-1996, revised
%           L. Mesin, 10-19-2016, revised to consider the translation of
%           half a sample
x=transpose(x(:));
x1=freshift(x,0.5); % shift of half a sample
vv=[x;x1];x=vv(1:end);
ntlag=2*ntlag;
n=length(x);
row=min(ntlag,n);
y=zeros(row,n/2);		% initialize the output matrix y
y(1,:)=x(1:2:end).*conj(x(1:2:end));	% zero lag instantaneous autocorrelation function
for i=1:row/2
    % positive lag instantaneous autocorrelation function
    v1=[x(i+1:n) zeros(1,i)];v2=conj([zeros(1,i) x(1:n-i)]);y(i+1,:)=v1(1:2:end).*v2(1:2:end);
end
i=fix(row/2);
v1=[x(i+1:n) zeros(1,i)];v2=conj([zeros(1,i) x(1:n-i)]);y(i+1,:)=v1(1:2:end).*v2(1:2:end);

% negative lag instantaneous autocorrelation function
y(end-row/2+2:end,:)=conj(flipud(y(2:row/2,:)));
% midtime lag instantaneous autocorrelation function
if i*2+1==row
    y(i+2,:)=conj(y(i+1,:));
    % if row is odd we also have a conj element of y(i+1,:)
end

end

function segt=freshift(seg,teta)
% Function to translate a signal in the frequency domain
SEG=fft(seg);
f=fftshift([-0.5:1/(length(seg)):0.5-1/(length(seg))]);
segt=ifft(SEG.*exp(1i*2*pi*teta*f));
end


