function [p,theta]=modelTPDmc(T,fitpar)
%function [p,theta]=modelTPDmc(T,fitpar)
%Model first order TPD curve as a sum of sites with coverages fitpar.thetas,
%fitpar should have following attributes:
%energies Ei=E0+i*dE,
%prefactor v,
%heating rate rate
%ml monolayer integral

% Returns:
% p - desorption rates array
% theta - coverages array

E0=fitpar.E0;
dE=fitpar.dE;
v=fitpar.v;
rate=fitpar.rate;
monolay=fitpar.ml;
thetas=fitpar.thetas;

k = 1.38e-23*6.24e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1

numpts=length(thetas);
E=linspace(E0,E0+dE*(numpts-1),numpts);
#Pre-calculate parameters
odepar.ek1=E/k; odepar.nuoa=v/rate;

if (rows(thetas)!=rows(T))
	thetas=thetas';
endif;

of=@(x,t) odemlde(x,t,odepar); #Solve the ODE
theta=lsode(of,thetas,T);

#Derive desorption rates from coverage values
numT=length(T);
p=zeros(numT,1);
for i=1:numT
  #For each temperature point find desorption flow as a sum of flows defined by remaining coverages
  p(i)=(-sum(odemlde(theta(i,:)',T(i),odepar)));
endfor

#Switch from 1/K to 1/sec.
p=p*rate*monolay;

endfunction;

#ODE for first order desorption
function dthetadT = odemlde(thetavec,T,param)
	nuovera=param.nuoa;
  Eoverk=param.ek1;
  dthetadT = (thetavec>eps)'.*(-nuovera).*thetavec'.*exp(-Eoverk./(T+eps));
endfunction;