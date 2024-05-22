% test the DCW estimate
close all;clear all
sigma=...  % choose the parameter of the Choi-Williams kernel
fc=256;		% sampling rate
T=1;		% window length
t=0:1/fc:T-1/fc;	% x axis
% Choose the frequency of two sine waves
fsin1=...;
fsin2=...;


% write different test signals (then uncomment a specific signal to make a test) 
% 1. sum of sinusoids 
A1=...;A2=...; % amplitudes of the sinusoids
y=A1*sin(2*pi*fsin1*t)+A2*sin(2*pi*fsin2*t);

% 2. sum of impulses (simulated as Gaussian functions with small std) 
alpha=.001;              % std of Gaussians approximating impulses 
% choose the time instants in which the two impulses are placed
tau1=...;
tau2=...;
y=...; % sum of 2 Gaussian functions centred in tau1 and tau2 and with std = alpha

% 3. chirp (simulate it using the function chirp)
y=chirp(t,...);
% Analytic signal: removing the contributions with negative frequency, less
% cross-terms are obtained
x=hilbert(y);					% analytic signal
% number of time lags to define tau
ntlag=...;	
% ambiguity function
aa=daf(...);
% definition of tau
tau=-ntlag/fc:1/fc:(ntlag-1)/fc;
% definition of theta
df=1/T;theta=-fc/2:df:fc/2-df;
% absolute value of the ambiguity function with the independent variables
% (tau,theta) reordered
af=fftshift(abs(aa));
% discrete Choi-Williams transformation
w=dcw(x,fc,ntlag,sigma,T);
% Wigner-Ville distribution
w1=wig(x,fc,ntlag);
% Frequency axis
[roww, colw]=size(w);
upper=fc/2;	% upper frequency
if 2*fix(roww/2)==roww; upper=upper-upper/roww; end
f=0:upper/(roww-1):upper;			% frequency axis

% Figure showing the results
figure,subplot(2,2,1);plot(t,y,'k')
subplot(2,2,2);contour(theta,tau,af,10,'k'),xlabel('Frequency lag \theta (Hz)'),ylabel('Time lag \tau (s)');title('Ambiguity function')
subplot(2,2,3);contour(t,f,real(w1)/max(max(real(w1))),.2:.2:1,'k'),xlabel('Time (s)'),ylabel('Frequency (Hz)');title('Wigner-Ville TFR')
subplot(2,2,4);contour(t,f,real(w)/max(max(real(w))),.2:.2:1,'k'),xlabel('Time (s)'),ylabel('Frequency (Hz)');title('Choi-Williams TFR')
