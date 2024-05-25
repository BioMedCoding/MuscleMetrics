function cv = beamopt3(Segna,posref,start,dint,fsamp);


 cpter=1;
 t=start;
 teta=10;
 trial=0; 

 while (abs(teta-t)>=5e-5) & trial < 30
  trial=trial+1;
  teta=t;     
  [de1, de2]=derivbeam(Segna,posref,teta);
  
  % Newton's criteria
  if (de2>0) 
    u=-de1/de2;
    if (abs(u)>0.5)
      u=-0.5*abs(de1)/de1;
    end
  else
    u=-0.5*abs(de1)/de1;
  end
  %err1(cpter,1)=de1;
  %err2(cpter,1)=de2;
  %u_v(cpter,1)=u;
  %cpter=cpter+1;
  
  %u=-de1/de2;
   
  t=teta+u;	  % result
 end

cv=dint/(teta/fsamp);