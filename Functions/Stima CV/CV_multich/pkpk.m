function cv=pkpk(Segna,dint,fsamp)

% Stima di CV studiando la distanza tra i picchi

for i = 1:1:size(Segna,1)
    [a b]=max(abs(Segna(i,:)));
    % parabola interpolante il segnale attorno al massimo
    A=[(b-1)^2 (b-1) 1;(b)^2 (b) 1; (b+1)^2 (b+1) 1];
    p=inv(A)*[Segna(i,b-1); Segna(i,b); Segna(i,b+1)]; % calcolo la parabola che passa per i 3 punti attorno al massimo 
    mx(i)=-p(2)/(2*p(1))/fsamp;
end

for i = 1:1:size(Segna,1)-1
estcv(i)=dint/(mx(i+1)-mx(i));    
end

cv=mean(estcv);