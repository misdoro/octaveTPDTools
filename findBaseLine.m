function ret=findBaseLine(tpd,debug=0)
di=[0;diff(tpd.i)];
dis=supsmu(tpd.T,di,'spa',0.01);

#Find the TPD half-maximum
maxi=max(find(tpd.i>max(tpd.i)/2));

#Find the baseline head end: min(i)
[mini,minidx]=min(tpd.i(1:maxi));
#Use minidx or at least first 10 points for the baseline head
headi=max(minidx,10);

#Find a di>0 point after the TPD half-maximum, this will be the baseline tail start
taili=[];
tpdlen=length(tpd.i);
dthresh=1e-11;
while(length(taili)==0 || taili > tpdlen-10)
  if (dthresh>eps)
    dthresh/=2;
  elseif (dthresh>0)
    dthresh=-eps;
  else
    dthresh*=2;
  endif
  posdis=find(dis>dthresh);
  taili=posdis(min(find(posdis>maxi+10)));
  if (debug)
    dthresh
    taili
    tpdlen
  endif
endwhile
#Works poorly on noisy TPDs

bl.T=cat(1,tpd.T(1:headi),tpd.T(taili:end));
bl.i=cat(1,tpd.i(1:headi),tpd.i(taili:end));

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