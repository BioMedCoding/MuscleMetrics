function x= ICWTW_LM(cwt, yscale,f,t)
%  Inverse Continuous Wavelet Transform
% Author: L. Mesin, 13/11/2016
%  Usage
%    x = ICWTW_LM(cwt, yscale,f,t)
%  Inputs
%    cwt      CWT
%    yscale   vector of scales
%    f      mother wavelet given in terms of a handle to a function
%    t        vector which defines the time axis of the signal
%  Outputs
%    x    signal
%
% Example
% global sd;Nlevels=1000;fsamp=128;sd=16/fsamp/10;N=4096;t=[1:N]/fsamp;
% sig=wfbm(.4,N);sig=sig-mean(sig);ft=.15;Wp=2*ft/fsamp;Rp=1;[b,a] = cheby1(5,Rp,Wp,'high');
% sig=filtfilt(b,a,sig);
% f=@(t) (1-(t-t(end/2+1)).^2/sd^2).*exp(-(t-t(end/2+1)).^2/2/sd^2);% Mexican hat
% [pxx_f,f_f] = pwelch(f(t),rectwin(length(sig)),[],[],fsamp);f_0=mean(f_f.*pxx_f)/mean(pxx_f);
% [pxx,fr] = pwelch(sig,rectwin(length(sig)),[],[],fsamp);
% vv=cumsum(pxx)/sum(pxx); fmin=max([fr(2) fr(max(find(vv<.01)))]);fmax=fr(min(find(vv>.99)));
% scale_min=f_0/fmax;scale_max=f_0/fmin;
% yscale=logspace(log10(scale_min),log10(scale_max),Nlevels);
% %yscale=linspace(scale_min,scale_max,Nlevels);
% [scalo, cwt] = CWTW_LM(sig,yscale,f,t);
% x = ICWTW_LM(cwt, yscale,f,t);err=std(x(:)-sig(:))/std(sig);
% subplot(2,1,1);plot(sig,'b');hold on;plot(x,'r');title('Input data and reconstruction');legend('input','after WT inversion');
% subplot(2,1,2);plot(x(:)-sig(:)); title(['Error (mean error = ' num2str(err*1e2,'%.1f') '%)'])
% 
dt=t(2)-t(1);
[n,nscale]= size(cwt);
% differential of the scale (used to compute the integral in the inversion formula)
vet=[2*yscale(1)-yscale(2) yscale 2*yscale(end)-yscale(end-1)];ds=(diff(vet(1:end-1))+diff(vet(2:end)))/2;
% Inizialization
x = zeros(n,1);
% Compute the amplitude scaling
fsamp=1/dt;t=t-t(end/2+1);df=1/max(t);Et=sum(abs(f(t)).^2)*dt;
xi =[-n/2:n/2-1]'*fsamp/n;dxi=xi(2)-xi(1);
window=fftshift(fft(f(t)));Ef=sum(abs(window).^2)*df;window=window(:)/sqrt(Ef)*sqrt(Et);
Cpsi=2*sum(abs(window(1:n/2)).^2./abs(xi(1:n/2)))*dxi;
for jo = 1:nscale
    % Current scale
    qscale = yscale(jo);
    % Scaled window
    tt=t/qscale;window_t=f(tt)/sqrt(qscale);window_t=window_t(:);    
    % Cross correlation integral on the delay
    w=conv(cwt(:,jo),window_t,'same')*dt;
    % Integral in the scale variable
    x = x+w/qscale^2*ds(jo)/Cpsi;
end
x=real(x);



