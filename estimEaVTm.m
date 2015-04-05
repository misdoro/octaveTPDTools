function Ea=estimEaVTm(v,Tm,rate)
Ea=0.01;
while (getV(Tm,Ea,rate)<v)
  Ea*=1.05;
endwhile;
endfunction;

function v=getV(Tm,Ea,rate)
k = 1.38e-23*6.24e18;
v=(rate*Ea.*exp(Ea/(k*Tm)))./(k*Tm.*Tm);
endfunction