#ODE for multilayer desorption, zero order, then first order
function dthetadt = odemlde(T,theta,param)
	nuovera=param.nuoa;
	Ek1=param.ek1;
	Eoverk=param.ek1;
	myzero =param.myzero;
	if (theta > 1)
		dthetadt = -nuovera.*exp(-Eoverk./(T+eps));
	elseif (theta > myzero)
		dthetadt = -nuovera.*theta.*exp(-Eoverk./(T+eps));
	else
		dthetadt = 0;
	end
endfunction;