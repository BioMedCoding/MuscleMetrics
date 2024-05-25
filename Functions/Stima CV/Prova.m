

clear;
close all
fs=2048;
t=0:1/fs:1;
ss=0.01;
sig1=diff(exp(-(t-.4).^2/2/ss^2));
del1=.001;
sig2=diff(exp(-(t-.4-del1).^2/2/ss^2));
sig3=diff(exp(-(t-.4-2*del1).^2/2/ss^2));
t=1/fs:1/fs:1;

ied=0.005;
sig=[sig1;sig2;sig3];
sig=sig/std(sig(:));
sig=sig+.1*randn(size(sig));

plot(t,sig(1,:));
hold on;
plot(t,sig(2,:)+1);
plot(t,sig(3,:)+2)
addpath('.\CV_multich')

cv_est = mle_CV_est(sig,ied,fs);