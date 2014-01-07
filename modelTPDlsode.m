function [ T, theta, p ] = modelTPDlsode(T,theta0,v,E,rate=3/60)
% Model TPD curve
%
% function [ T, theta, p ] = modelTPDlsode(T,theta0,v,E,rate=3/60)
%
% T - temperatures needed, array of points
% theta0 - initial coverage
% v - frequency (desorption prefactor)
% E - adsorption energy in electron-volts
% rate - desorption rate, defaults to 3K/minute
%
% Returns:
% T - the temperatures array
% theta - coverages array
% p - desorption rates array


k = 1.38e-23*6.24e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1
#Pre-calculate parameters
odepar.ek1=E/k; odepar.nuoa=v/rate; odepar.myzero=eps;

of=@(x,t) odemlde(x,t,odepar); #Solve the ODE
theta=lsode(of,theta0,T);

if (rows(theta)!=rows(T))
	theta=theta';
endif;

#Derive desorption rates from the coverage values
p=-arrayfun(@odemlde,theta,T,odepar);

endfunction;

#ODE for multilayer desorption, zero order, then first order
function dthetadT = odemlde(theta,T,param)
	nuovera=param.nuoa;Eoverk=param.ek1;myzero=param.myzero;
	if (theta > 1)
		dthetadT = -nuovera.*exp(-Eoverk./(T+eps));
	elseif (theta > myzero)
		dthetadT = -nuovera.*theta.*exp(-Eoverk./(T+eps));
	else
		dthetadT = 0;
	end
endfunction;