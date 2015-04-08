#Estimate number of point to cover the used temperature range
function np=estimNP(mytpd,fit)
#Find np such that max(mytpd.T)-1=Tmax(E0+np*dE)
emax=estimEaVTm(fit.v,max(mytpd.T)-1,fit.rate);
np=round((emax-fit.E0)/fit.dE);
#input("break");
endfunction
