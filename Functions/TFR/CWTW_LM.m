function [scalo, cwt] = CWTW_LM(x,yscale,f,t)
% Continuous Wavelet Transform
% Author: L. Mesin, 13/11/2016
%  Usage
%    [scalo, cwt] = CWTW(x,yscale,wavelet,t)
%  Inputs
%    x        signal, dyadic length n=2^J, real-valued
%    yscale   vector of scales
%    f        mother wavelet given in terms of a handle to a function
%    t        vector which defines the time axis of the signal
%  Outputs
%    scalo    Scalogram. Matrix n by nscale where
%    cwt      CWT
% Example
% global sd;Nlevels=1000;fsamp=128;sd=16/fsamp/10;N=4096;t=[1:N]/fsamp;
% sig=wfbm(.4,N);sig=sig-mean(sig);ft=.15;Wp=2*ft/fsamp;Rp=1;[b,a] = cheby1(5,Rp,Wp,'high');
% sig=filtfilt(b,a,sig);
% [pxx,fr] = pwelch(sig,rectwin(length(sig)),[],[],fsamp);
% f=@(t) (1-(t-t(end/2+1)).^2/sd^2).*exp(-(t-t(end/2+1)).^2/2/sd^2);% Mexican hat
% [pxx_f,f_f] = pwelch(f(t),rectwin(length(sig)),[],[],fsamp);f_0=mean(f_f.*pxx_f)/mean(pxx_f);
% vv=cumsum(pxx)/sum(pxx); fmin=max([fr(2) fr(max(find(vv<.01)))]);fmax=fr(min(find(vv>.99)));
% scale_min=f_0/fmax;scale_max=f_0/fmin;
% yscale=logspace(log10(scale_min),log10(scale_max),Nlevels);
% %yscale=linspace(scale_min,scale_max,Nlevels);
% [scalo, cwt] = CWTW_LM(sig,yscale,f,t);
% figure;subplot(2,1,1);plot(t,sig);title('Data');xlabel('Time (s)')
% subplot(2,1,2);image(t,yscale,(flipud(scalo')-min(scalo(:)))/range(scalo(:))*256);
% colormap(pink(255)); xlabel('Time (s)');ylabel('Scales')

dt=t(2)-t(1);
n = length(x);
% Inizialization
nscale=length(yscale);
cwt = zeros(n,nscale);scalo=cwt;
kscale  = 1;
for jo = 1:nscale
    % Compute the current scale
    qscale = yscale(jo);
    % Compute the scaled window
    tt=t/qscale;window_t=f(tt)/sqrt(qscale);window_t=window_t';window_t=window_t(:);
    % Flip of the signal to get a cross - correlation
    % using the conv command
    % (the first sample is discarded in order that 
    % the window remains centred at end/2+1).
    window_t=flipud(window_t(3:end));
    % Convolution integral
    w=conv(x,window_t,'same')*dt;
    % Scalogram
    scalo(1:n,kscale) = (abs(w)).^2;
    % CWT
    cwt(1:n,kscale) = w;
    kscale = kscale+1;
end

