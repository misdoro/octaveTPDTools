function [ Tode, theta, p ] = modelTPDlsode(T,theta0,v,E,rate=3/60)
% Model TPD curve
%
% T - temperatures needed, independent variable
% N - normalization factor, area under p(T)
% v - frequency
% E - energy level in electron volts


k = 1.38e-23*6.24e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1
odepar.ek1=E/k;
odepar.nuoa=v/rate;
odepar.myzero=10*eps;

global odparam=odepar;


tic();
theta=lsode(@odemlde,theta0,T);
printf("lsode time:")
toc()


#if (isfield(ret,"stats"))
#	ret.stats
#endif
if (rows(T)>1)
	Tode=T;
else
	Tode=T';
endif;
#theta=ret.y;
p=-arrayfun(@odemlde,theta,Tode);

endfunction;

#ODE for multilayer desorption, zero order, then first order
function dthetadT = odemlde(theta,T)
	global odparam;
	param=odparam;
	nuovera=param.nuoa;
	Ek1=param.ek1;
	Eoverk=param.ek1;
	myzero =param.myzero;
	if (theta > 1)
		dthetadT = -nuovera.*exp(-Eoverk./(T+eps));
	elseif (theta > myzero)
		dthetadT = -nuovera.*theta.*exp(-Eoverk./(T+eps));
	else
		dthetadT = 0;
	end
endfunction;