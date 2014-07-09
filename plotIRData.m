function result=plotIRData(mytpd,param,result,press,dose,varargin)
%  IRDAT=loadIRData();
%  figure(1);
%  clf();
%  hold on;
%  colors=linspace(1,4,IRDAT.count);
%  for i=1:IRDAT.count;
%    data=IRDAT.data{i};
%    plot(data.wl,data.abs,'color',getLineColor(colors(i),1,2));
%  endfor;
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

for i=1:irdat.count
  irtime=mktime(irdat.data{i}.time);
  if (irtime>=tpd.mintime && irtime<=tpd.maxtime)
    count++;
    irplot{count}=irdat.data{i};
    irtimes(count)=irtime;
  endif
endfor

irpoints={};
if (count>0)
  for i=1:count
    color=getLineColor(4*i/count,1,2);
    figure(param.fig.IR);
    plot(irplot{i}.wl,irplot{i}.abs,'color',color);
    if (isfield(param.fig,"disp"))
      figure(param.fig.disp);
      tpdtime=mktime(mytpd.time);
      irt=irtimes(i)-tpdtime;
      idx=min(find(mytpd.t>irt));
      plot(mytpd.T(idx),mytpd.i(idx),"x",'linewidth',2,"color",color);
      
    endif
    
  endfor
endif
result=irpoints;

endfunction;