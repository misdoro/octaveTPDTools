function result=plotIRData(mytpd,param,result,press,dose,varargin)
% Some initial tests
if (length(varargin)<1)
  printf("Need an IRDAT argument\n");
  return;
endif

if (~isfield(mytpd,'time'))
  printf("TPD has no start time field\n");
  return;
endif

irdat=varargin{1}{1};
if (~isstruct(irdat) || ~isfield(irdat,'count') || ~irdat.count>0)
  printf("The supplied argument is not the IR data structure or has no data\n");
endif


tpd.mintime=mktime(mytpd.time)+min(mytpd.t);
tpd.maxtime=mktime(mytpd.time)+max(mytpd.t);
count=0;
irplot={};

if isfield(param,"IRTemp")
  IRTemp=sort(param.IRTemp);
else
  IRTemp=[];
endif

for i=1:irdat.count
  irtime=mktime(irdat.data{i}.time);
  if (irtime>=tpd.mintime && irtime<=tpd.maxtime)
    if ((length(IRTemp) && checkIRTemp(IRTemp,mytpd,irtime))||~length(IRTemp))
      count++;
      irplot{count}=irdat.data{i};
      irtimes(count)=irtime;
    endif
  endif
endfor

irpoints={};
if (count>0)
  for i=1:count
    color=getLineColor(4*i/count,1,2);
    figure(getFigIndex("IR"));
    if (isfield(param,"IRShift"))
      if length(param.IRShift)==1
        irplot{i}.abs+=param.IRShift*i;
      else
        irplot{i}.abs+=param.IRShift(1)+param.IRShift(2)*i;
      endif
    endif
    plot(irplot{i}.wl,irplot{i}.abs,'color',color);
    if (getFigIndex("disp",0))
      figure(getFigIndex("disp"));
      tpdtime=mktime(mytpd.time);
      irt=irtimes(i)-tpdtime;
      idx=min(find(mytpd.t>irt));
      plot(mytpd.T(idx),mytpd.i(idx),"x",'linewidth',2,"color",color);
      
    endif
    
  endfor
endif
result=irpoints;
endfunction;

function IRT = getIrTemp(mytpd,irtime)
  timeidx=min(find(mytpd.t>(irtime-mktime(mytpd.time))));
  IRT=mytpd.T(timeidx);
endfunction;

function ret=checkIRTemp(IRTemp,mytpd,irtime)
  ret=find(abs(IRTemp-getIrTemp(mytpd,irtime))<2);
endfunction