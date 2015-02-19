function ret=saveFig(label,fileprefix,saveAscii=0)
  %A function to save a figure as image and its plot data as ascii
  hdl=getFigIndex(label,0);
  if (hdl)
    %Save figure as image
    print(hdl,fileprefix,"-dpng","-r300");
    printf("Saved %s.png\n",fileprefix);
    %Extract and save figure data
    if (saveAscii)
      saveFigData(hdl,fileprefix);
    endif;
  endif
endfunction;

function saveFigData(fighdl,fileprefix);
  graphs=get(fighdl,'children');
  plots=get(graphs(length(graphs)),'children');
  idx=0;
  for i=1:length(plots);
    refi=plots(i);
    objtype=get(refi,'type');
    if (strcmp(objtype,'line'))
      saveCurveData(refi,fileprefix,++idx);
    endif;
  endfor
endfunction;

function saveCurveData(plotref,fprefix,index)
  xdata=get(plotref,'xdata')';
  ydata=get(plotref,'ydata')';
  fname=sprintf("%s.%02d.asc",fprefix,index);
  plotdata=[xdata,ydata];
  save("-ascii",fname,'plotdata');
endfunction;