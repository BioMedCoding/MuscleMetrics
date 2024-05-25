function [dist, scale, b]=sig_dis(signal1,signal2,fsamp,dista)

% Final version by Dario Farina 
% 1/3/2000

long=length(signal1);
signal1=abs(signal1);
signal2=abs(signal2);

    X1=cumsum(signal1);            % La function de repartition

    t=linspace(1/1024,1,long);
         X2=X1./max(X1);              % Normalisation
         for JJ=1:long-1,
            if X2(JJ)>=X2(JJ+1)
            X2(JJ+1)=X2(JJ)+0.0001;
            end
         end

	Y1=cumsum(signal2);   % F. R.

   Y2=Y1./max(Y1);             % Normalisation
         for JJ=1:long-1,
            if Y2(JJ)>=Y2(JJ+1)
            Y2(JJ+1)=Y2(JJ)+0.0001;
            end
         end
         
         tt=interp1(X2,t,[0.01:0.01:0.99]);
         
         ttap=interp1(Y2,t,[0.01:0.01:0.99]);
         
         a1=find(~isnan(tt));
         
         a2=find(~isnan(ttap));
         e=1;
         for i=1:length(a1)
            if ~isempty(find(a1(i)==a2))
               atot(e)=a1(i);
               e=e+1;
            end;
         end;
         
         ttap=ttap(atot);
         tt=tt(atot);
         
         [P S]=polyfit(ttap,tt,1);       
         %close all
         %plot(ttap,tt,'k')
         %hold on
         %plot(ttap,ttap+P(2),'r')
         %plot(ttap,P(1)*ttap+P(2),'g')         
         %pause
  
  b=P(2)*length(signal1)/fsamp;
  
  dist=0;
  
  for i=1:length(tt)
     if dista == 'dist_scale'
        dist=dist+abs(tt(i)-(ttap(i)+P(2)))/sqrt(2);
     else
        dist=dist+abs(tt(i)-(ttap(i)*P(1)+P(2)))/sqrt(2);
     end;        
  end;
  
  dist=dist/length(tt);
  
  scale = P(1);
     
     
 

