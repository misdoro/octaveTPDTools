function ret=removeSpikes(x,y)
  ysm=supsmu(x,y,'spa',0.005);
  dy=y-ysm;
  stddy=std(dy)*10;
  sp=find(abs(dy)>stddy);
  for i=1:length(sp);
    spi=sp(i);
    prevy=y(spi);
    newy=mean([y(spi-1),y(spi+1)]);
    y(spi)=newy;
    printf("replaced a spike at %d of %.3e with %.3e\n",spi,prevy,newy);
  end
  ret=y;
endfunction