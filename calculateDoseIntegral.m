function [doseint,maxidx]=calculateDoseIntegral(dosedata,debug=0)
%Calculate dose integral: extract the baseline, finish integration at 1/2(max-end)
%Returns the dose integral and the final point of the integration
	
	if (!isstruct(dosedata) ||
		!isfield(dosedata,"t")||
		!isfield(dosedata,"i")||
		length(dosedata.i)==0 ||
		length(dosedata.t)==0);
		doseint=0;
		maxidx=0;
		return;
	endif;

	dose.mini=min(dosedata.i);
	[dose.maxi,dose.maxii]=max(dosedata.i);
  dose.decay=dosedata.i(dose.maxii:end);
  dose.miniend=min(dose.decay);
	dose.fini=dosedata.i(end);
	dose.iend=(dose.maxi+2*dose.miniend)/3;
  #dose.iend=dose.maxi;
	dose.intend=dose.maxii+min(min(find(dose.decay<=dose.iend)),length(dosedata.i));
	dose.isub=dosedata.i(1:dose.intend)-dose.mini;
	dose.tsub=dosedata.t(1:dose.intend);
	doseint=trapz(dose.tsub,dose.isub);
	maxidx=dose.intend;
  
  if(debug)
    figure(99)
    clf();
    hold on;
    plot(dosedata.t,dosedata.i);
    plot(dosedata.t(dose.intend),dosedata.i(dose.intend),'x','color','red');
	  input("OK?")
  endif;
	
  endfunction;
