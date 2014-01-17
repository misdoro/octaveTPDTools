#Fit temperature slope with linear function, interpolate points
function fittemp = fitTempSlope(timepoints,tempdata,Trange)
	data.t=timepoints;
	data.T=tempdata;
	dcut=cutTemp(data,Trange);
	
	fittemp.poly=polyfit(dcut.t,dcut.T,1);
	fittemp.points=polyval(fittemp.poly,timepoints);
endfunction