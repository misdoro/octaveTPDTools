function ret=removeSpikes(x,y)
  ysm=supsmu(x,y,'spa',0.005);
  dy=y-ysm;
  stddy=std(dy)*10;
  sp=find(abs(dy)>stddy);
  sze=length(x);
  for i=1:length(sp);
    spi=sp(i);
    if (spi<sze && spi>1)%Proceed only if spike is detected before the last element and after the first.
      prevy=y(spi);
      newy=mean([y(spi-1),y(spi+1)]);
      y(spi)=newy;
      printf("replaced a spike at %d of %.3e with %.3e\n",spi,prevy,newy);
    endif
  end
  ret=y;
endfunction