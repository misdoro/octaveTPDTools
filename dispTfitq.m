function result=dispTfitq(mytpd,param,result,press,dose);

figure(getFigIndex("Tfitq"));
plot(press.T,(press.To-press.T),"linewidth",1,"color",mytpd.color);
figure(getFigIndex("Tramp"));
minT=min(press.To);
maxT=press.To(end);
rampT=min(max(press.T,minT),maxT);
plot(press.t,press.To,press.t,rampT);
endfunction;