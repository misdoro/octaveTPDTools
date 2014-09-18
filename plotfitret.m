function ret=plotfitret(filename,linecolor='blue');
load(filename);
ret.v=fitret(:,1);
ret.Ea=fitret(:,2);
ret.cov=fitret(:,3);
ret.ssq=fitret(:,4);

figure(1);
hold on;
semilogx(ret.v,ret.Ea,'color',linecolor);
xlabel("v");
ylabel("Ea");

figure(2);
hold on;
semilogx(ret.v,ret.ssq/min(ret.ssq),'color',linecolor);
xlabel("v");
ylabel("fit rms");

figure(3);
hold on;
semilogx(ret.v,ret.cov,'color',linecolor);
xlabel("v");
ylabel("fit coverage");
endfunction;