function blval=baselineFindConst(mytpd); 
%Function to find a constant baseline in TPD record
  [mini,minidx]=min(mytpd.i);
  blst=max(minidx-20,1);
  blend=min(minidx,length(mytpd.i));
  blval=mean(mytpd.i(blst:blend));
endfunction;