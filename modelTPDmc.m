function [p,theta]=modelTPDmc(T,thetavec,par)
%function [theta,p]=modelTPDmc(T,thetavec,par)
%Model TPD curve as a sum of sites with coverages thetavec,
%par has par.E0,par.dE,par.v,par.rate attributes
%energies Ei=E0+i*dE,
%prefactor v,
%heating rate rate

% Returns:
% p - desorption rates array
% theta - coverages array

E0=par.E0;
dE=par.dE;
v=par.v;
rate=par.rate;
monolay=par.monolay;

k = 1.38e-23*6.24e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1

numpts=length(thetavec);
E=linspace(E0,E0+dE*(numpts-1),numpts);
#Pre-calculate parameters
odepar.ek1=E/k; odepar.nuoa=v/rate;

if (rows(thetavec)!=rows(T))
	thetavec=thetavec';
endif;

of=@(x,t) odemlde(x,t,odepar); #Solve the ODE
theta=lsode(of,thetavec,T);

#Derive desorption rates from the coverage values
numT=length(T);
p=zeros(numT,1);

#This one is slow...
for i=1:numT
p(i)=(-sum(odemlde(theta(i,:)',T(i),odepar)));
endfor

p=p*rate*par.monolay;

endfunction;


#ODE for first order desorption
function dthetadT = odemlde(thetavec,T,param)
	nuovera=param.nuoa;
  Eoverk=param.ek1;
  dthetadT = (thetavec>eps)'.*(-nuovera).*thetavec'.*exp(-Eoverk./(T+eps));
endfunction;
function dthetadT = odemlde1(thetavec,T,param)
	nuovera=param.nuoa;
  Eoverk=param.ek1;
  Eoverk=Eoverk(param.i);
  dthetadT = (thetavec>eps).*(-nuovera).*thetavec.*exp(-Eoverk./(T+eps));
endfunction;