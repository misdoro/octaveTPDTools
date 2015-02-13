function [figidx]=getFigIndex(figlabel,create=1);
%A function to get a figure number from its text label
	persistent figindex=0;
  persistent figs=struct();
	if (isfield(figs,figlabel))
		figidx=getfield(figs,figlabel);
	elseif(create)
		figidx=++figindex;
		figs=setfield(figs,figlabel,figidx);
  else
    figidx=0;
	endif;
return;