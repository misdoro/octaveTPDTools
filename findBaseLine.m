function ret=findBaseLine(tpd,debug=0)
if (length(tpd.i)!=length(tpd.t))
  return
endif
di=[0;diff(tpd.i)];
dis=supsmu(tpd.i,di,'spa',0.01);

#Minimum number of points at each end
minpts=10;
#Find the TPD half-maximum
maxi=max(find(tpd.i>max(tpd.i)/2));

#Find the baseline head end: min(i)
[mini,minidx]=min(tpd.i(1:maxi));
#Use minidx or at least first 10 points for the baseline head
headi=max(minidx,minpts);

tpdlen=length(tpd.i);
taili=tpdlen;
if (maxi+2*minpts>tpdlen)
  taili=tpdlen-minpts;#If tpd half-maximum is way too close or 
  # behind the end, take the last minpts points
else
  #Find a di>0 point after the TPD half-maximum, 
  #this will be the baseline tail start
  dthresh=1e-11;
  while(length(taili)==0 || taili > tpdlen-minpts)
    if (dthresh>eps)
      dthresh/=2;
    elseif (dthresh>0)
      dthresh=-eps;
    else
      dthresh*=2;
    endif
    posdis=find(dis>dthresh);
    taili=posdis(min(find(posdis>maxi+minpts)));
    if (debug)
      dthresh
      taili
      tpdlen
    endif
  endwhile
endif
#Works poorly on noisy TPDs


bl.t=cat(1,tpd.t(1:headi),tpd.t(taili:end));
bl.i=cat(1,tpd.i(1:headi),tpd.i(taili:end));

[poly,s]=polyfit(bl.t,bl.i,3);
tpd.bl=polyval(poly,tpd.t);
tpd.ipur=max(tpd.i-tpd.bl,0);

ret=tpd;
if (debug)
  figure(1)
  plot(tpd.t,dis,tpd.t,di)
  figure(2)
  clf()
  hold on
  plot(tpd.t,tpd.i)
  plot(tpd.t,polyval(poly,tpd.t),'color',"green")
  plot(bl.t,bl.i,'color',"red")
  plot(tpd.t,tpd.ipur);
  input("OK?");
endif