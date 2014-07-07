function result=dispTfitq(mytpd,param,result,press,dose);

#plot(press.t,press.To,"linewidth",2,"color",mytpd.color);
#plot(press.t,press.T,"linewidth",1,"color",mytpd.color,'linestyle',':');
plot(press.T,(press.To-press.T),"linewidth",1,"color",mytpd.color);

endfunction;