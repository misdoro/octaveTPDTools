function Ea=estimEaVTm(v,Tm,rate)
Ea=0.01;
while ((dv=(getV(Tm,Ea,rate)-v))<0)
  Ea*=1.1;
endwhile;
while ((dv=(getV(Tm,Ea,rate)-v))>0)
  Ea*=0.999;
endwhile;
endfunction;

function v=getV(Tm,Ea,rate)
k = 1.38e-23*6.24e18;
v=(rate*Ea.*exp(Ea/(k*Tm)))./(k*Tm.*Tm);
endfunction