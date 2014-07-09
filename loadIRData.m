function ret=loadIRData(path='')
fid=fopen(strcat(path,'time.txt'));
line='';
ircount=0;
if (fid>0)
  while((line=fgetl(fid))~=-1);
    irfile=line;
    tline=strcat(fgetl(fid),irfile);
    irtime=strptime(tline,"%H:%M:%SIR_%Y%m%d");
    asctime(irtime);
    ircount++;
    tab=dlmread(strcat(path,irfile));
    datpt.time=irtime;
    datpt.wl=tab(:,1);
    datpt.abs=tab(:,2);
    irdat{ircount}=datpt;
  endwhile
else
ircount=0;
irdat=[];
endif
ret.count=ircount;
ret.data=irdat;
endfunction