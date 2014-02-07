function result=calculateDoseIntegral(dosedata)
%Calculate dose integral: extract the baseline, finish integration at 1/2(max-end)
	dose.mini=min(dosedata.i);
	dose.maxi=max(dosedata.i);
	dose.fini=dosedata.i(end);
	dose.iend=(dose.maxi+dose.fini)/2;
	dose.intend=min(max(find(dosedata.i>dose.iend)),length(dosedata.i));
	dose.isub=dosedata.i(1:dose.intend)-dose.mini;
	dose.tsub=dosedata.t(1:dose.intend);
	#figure(100);
	#plot(dose.tsub,dose.isub);
	#drawnow();
	#input("OK?");
	result=trapz(dose.tsub,dose.isub);

endfunction;
