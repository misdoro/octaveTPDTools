function [ Tode, theta, p ] = modelTPDcd(T,theta0,v,E,rate=3/60,alpha)
% Model TPD curve, coverage-dependant
%
% T - temperatures needed, independent variable
% N - normalization factor, area under p(T)
% v - frequency
% E - energy level in electron volts


k = 1.38e-23*6.24e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1
odepar.ek1=E/k;
odepar.nuoa=v/rate;
odepar.myzero=10*eps;
odepar.alphaok=alpha/k;


vopt=odeset("NormControl","on","Stats","off");

tic();
ret=ode45(@odemlde,T,theta0,vopt,odepar);
toc()

if (isfield(ret,"stats"))
	ret.stats
endif
Tode=ret.x;
theta=ret.y;
odepar
p=-arrayfun(@odemlde,Tode,theta,odepar);

endfunction;

#ODE for multilayer desorption, zero order, then first order
function dthetadt = odemlde(T,theta,param)
	nuovera=param.nuoa;
	Ek1=param.ek1;
	Eoverk=param.ek1;
	myzero =param.myzero;
	alphaok=param.alphaok;
	if (theta > 1)
		dthetadt = -nuovera.*exp(-Eoverk./(T+eps));
	elseif (theta > myzero)
		dthetadt=-nuovera.*theta.*exp(-(Eoverk-alphaok*theta)./(T+eps));
	else
		dthetadt = 0;
	end
endfunction;

function ratio=getRatio(theta,param)
	
endfunction;
