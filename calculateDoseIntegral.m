function [doseint,maxidx]=calculateDoseIntegral(dosedata)
%Calculate dose integral: extract the baseline, finish integration at 1/2(max-end)
%Returns the dose integral and the final point of the integration
	
	if (!isstruct(dosedata) ||\
		!isfield(dosedata,"t")||\
		!isfield(dosedata,"i")||\
		length(dosedata.i)==0||\
		length(dosedata.t)==0);
		doseint=0;
		maxidx=0;
		return;
	endif;
	
	dose.mini=min(dosedata.i);
	dose.maxi=max(dosedata.i);
	dose.fini=dosedata.i(end);
	dose.iend=(dose.maxi+dose.fini)/2;
	dose.intend=min(max(find(dosedata.i>dose.iend)),length(dosedata.i));
	dose.isub=dosedata.i(1:dose.intend)-dose.mini;
	dose.tsub=dosedata.t(1:dose.intend);
	doseint=trapz(dose.tsub,dose.isub);
	maxidx=dose.intend;
endfunction;
