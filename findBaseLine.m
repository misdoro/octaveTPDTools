function ret=findBaseLine(tpd,debug=0,dthresh)
di=[0;diff(tpd.i)];
dis=supsmu(tpd.T,di,'spa',0.01);

#Find a di>0 point after the TPD maximum, this will be the baseline start
[maxv,maxi]=max(tpd.i);
dislim=[0 1e-13 2e-13 5e-13 1e-12];
posdis=find(dis>dthresh);
bltailidx=posdis(min(find(posdis>maxi+10)));

bl.T=cat(1,tpd.T(1:10),tpd.T(bltailidx:end));
bl.i=cat(1,tpd.i(1:10),tpd.i(bltailidx:end));

[poly,s]=polyfit(bl.T,bl.i,3);
tpd.bl=polyval(poly,tpd.T);
tpd.ipur=max(tpd.i-tpd.bl,0);

ret=tpd;
if (debug)
  figure(1)
  plot(tpd.T,dis,tpd.T,di)
  figure(2)
  clf()
  hold on
  plot(tpd.T,tpd.i)
  plot(tpd.T,polyval(poly,tpd.T),'color',"green")
  plot(bl.T,bl.i,'color',"red")
  plot(tpd.T,tpd.ipur);
endif