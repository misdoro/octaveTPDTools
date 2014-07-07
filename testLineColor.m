yv=linspace(0,10,100);
x=[0,1];
figure(2);
clf();
hold on;
for y=yv
  cl=getLineColor(y,1,2);
  plot(x,[y,y],'color',cl);
endfor