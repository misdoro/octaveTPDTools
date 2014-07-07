function ret=loadIRData()
fid=fopen('time.txt');
line='';
ircount=0;

while((line=fgetl(fid))~=-1);
  irfile=line;
  tline=strcat(fgetl(fid),irfile);
  irtime=strptime(tline,"%H:%M:%SIR_%Y%m%d");
  asctime(irtime);
  ircount++;
  tab=dlmread(irfile);
  datpt.time=irtime;
  datpt.wl=tab(:,1);
  datpt.abs=tab(:,2);
  irdat{ircount}=datpt;
endwhile
ret.count=ircount;
ret.data=irdat;
endfunction