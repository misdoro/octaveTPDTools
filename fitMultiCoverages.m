function result=fitMultiCoverages(mytpd,param,result);
  if (~isfield(result,"fitcov"))
    #Load fit data from file on first iteration
    [f.info, f.err, f.msg]=stat("fitcov.par");
    if (f.err>=0);
      load "fitcov.par";
      result.fitcov=fitcov;
      fitde=findFileMetadata(result.fitcov,mytpd.filename);
    else
      fitde.found=0;
    endif
  else
    fitde=findFileMetadata(result.fitcov,mytpd.filename);
  endif
      
  if (~fitde.found)
    if (isfield(result,"fitcov"))
      printf("Using parameters from previous fit iteration\n");
      fitpar=result.fitcov{mytpd.idx-1};
      fitpar.rate=mytpd.rate;
    else
      fitpar=myInitFitPar(mytpd,param)
    endif
    if (isfield(param,"fitpenalty"))
      fitpar.penalty=param.fitpenalty;
    endif
    if (isfield(param,"fitmaxiter"))
      maxiter=param.fitmaxiter;
    else
      maxiter=1500;
    endif
    if (isfield(param,"debug"))
      fod="iter"
    else
      fod="final"
    endif
    fitopts=optimset("Display",fod,"MaxIter",maxiter,"TolX",1e-5);
    fitde=fitPartCoverages(mytpd,fitpar,fitopts);
    fitde.filename=mytpd.filename;
  else
    printf("Plotting saved fit data\n");
  endif
  
  np=length(fitde.thetas);
  Epts=linspace(fitde.E0,fitde.E0+(np-1)*fitde.dE,np);
  
  printf("Total initial coverage %.2f ML\n",sum(fitde.thetas));
  [wmeanea,efwhm]=getemaxde(Epts,fitde.thetas);
  printf("Mean weighted Ea=%.3f eV, dE(FWHM)=%.3f eV\n",wmeanea,efwhm);
  fitrow=[sum(fitde.thetas),wmeanea,efwhm];
  result=retAppend(result,"fitemaxde",fitrow);
  
  figure(getFigIndex("covsites"));
  plot(Epts,fitde.thetas,"linewidth",2,"color",mytpd.color);
  
  p=modelTPDmc(mytpd.T,fitde);
  #Check if we are asked to convert the data in ML/K units
  if (index(param.tools,'N'))
    p=p./(param.monolayer.*(mytpd.rate));
  endif
  figure(getFigIndex("disp"))
  plot(mytpd.T,p,"color",mytpd.color,"linestyle","--");
  figure(getFigIndex("covfiterr"))
  plot(mytpd.T,10*(mytpd.fileidx-1)+100*(p-mytpd.i)/max(mytpd.i),"color",mytpd.color);
  drawnow();
  
  fcl=1;
  if (isfield(result,"fitcov"))
    if (~(fcl=findFileMetadata(result.fitcov,mytpd.filename).found))
      fcl=length(result.fitcov)+1;
    endif;
  endif;
  result.fitcov{fcl}=fitde;
  fitcov=result.fitcov;
  save("-text","fitcov.par","fitcov");
  
endfunction;

function ret=findFileMetadata(cellarr,filename)
  ret=struct("found",0);
  for i=1:length(cellarr)
    if (isfield(cellarr{i},"filename") && cellarr{i}.filename==filename)
      ret=cellarr{i};
      ret.found=i;
      break;
    endif
  endfor
    
endfunction

function fitpar=myInitFitPar(mytpd,param)
  fit=loadFitFile();
  fitpar=struct();
  if (isfield(fit,"fitdE"))
    printf("Using best fit prefactor %.1e from fit.par\n",fit.fitdE.v);
    fitpar=initFitParam(mytpd,param,fit.fitdE.v);
  else
    fitpar=initFitParam(mytpd,param);
  endif;
  
  offset=6;
  fitpar.E0-=offset*fitpar.dE;
  fitpar.thetas=0.01*ones(fitpar.np+offset,1);
  
  fitpar.rate=mytpd.rate;
endfunction

function [avgx] = wmean(x,weight,wpow=4)
  wx = x * weight.^wpow;
  avgx = sum (wx) / sum (weight.^wpow);
endfunction

function [emean, efwhm]=getemaxde(Epts,thetas)
  emean=wmean(Epts,thetas);
  eptsus=linspace(Epts(1),Epts(end),100);
  covus=interp1(Epts,thetas,eptsus);
  covmax=covus(min(find(eptsus>emean)));
  maxe=eptsus(max(find(covus>covmax/2)));
  mine=eptsus(min(find(covus>covmax/2)));
  efwhm=maxe-mine;
endfunction