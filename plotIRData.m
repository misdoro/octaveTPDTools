function plotIRData()
  IRDAT=loadIRData();
  figure(1);
  clf();
  hold on;
  colors=linspace(1,4,IRDAT.count);
  for i=1:IRDAT.count;
    data=IRDAT.data{i};
    plot(data.wl,data.abs,'color',getLineColor(colors(i),1,2));
  endfor;
endfunction;