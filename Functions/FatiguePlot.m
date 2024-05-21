function [rms,arv,mnf,mdf,cv]=FatiguePlot(sig,s1,s2,fs,IED)

%misure di ampiezza
rms=sqrt(mean(sig.*sig));   %rms
arv=sum(abs(sig))/length(sig);  %arv

%misure di frequenza
[Pxx,f] = pwelch(sig-mean(sig),tukeywin(length(sig),.1),0,[],fs);
mnf=sum(f.*Pxx)/sum(Pxx); %frequenza media MNF

cumulative_psd=cumsum(Pxx);
median_value= cumulative_psd(end)/2;
median_index = find(cumulative_psd >= median_value, 1);
mdf = f(median_index); %frequenza mediana MDF

%velocità di conduzione
[xc,del]=xcorr(s2,s1,10);
[mm,I]=max(xc);
start=del(I);
d=delay(real(fft(s2)),imag(fft(s2)),real(fft(s1)),imag(fft(s1)),start);
cv=IED/(d*1000)*fs;
end