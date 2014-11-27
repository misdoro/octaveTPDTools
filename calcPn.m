function ptotn=calcPn(tpd,parv,par)
	%TODO:To throw away, too messy!
	#Sum of tpds with variable parameters
	numsim=par.numsim;
	pars=parv;
	if (~isfield(par,"thetai") || length(par.thetai)~=numsim)
		pars(3)=parv(3)/numsim;
	else
		pars(3)=parv(3);
	endif;
	parm=repmat(pars,numsim,1);
  if (length(parv)==5)
    espr=parv(5);
    cspr=0;
  else
  	espr=par.Es;
	  cspr=par.covspr;
	endif
  if (isfield(par,"Ei") && length(par.Ei)==numsim);
		parm(:,4)=par.Ei;
	else
		parm(:,4)=linspace(pars(4)-espr,pars(4)+espr,numsim);
	endif;
	if (isfield(par,"thetai") && length(par.thetai)==numsim);
		parm(:,2)=par.thetai;
	else
		parm(:,2)=linspace(pars(2)*(1-cspr),pars(2)*(1+cspr),numsim);
	endif;
	
	
	if (par.parallel<2);
		ptotn=zeros(length(tpd.T),1);
	
		for i=1:numsim
			if (isfield(par,"ls")&&par.ls==1);
				ptotn+=calcP1(tpd,parm(i,:));
			else
				ptotn+=calcP1o(tpd,parm(i,:));
			endif
		endfor
	else
	
    pkg load parallel;
		for i=1:numsim
			para(i).a=parm(i,1);
			para(i).b=parm(i,2);
			para(i).c=parm(i,3);
			para(i).d=parm(i,4);
			para(i).T=tpd.T;
			para(i).rate=tpd.rate;
		endfor
		if (isfield(par,"ls")&&par.ls==1);
			ptotn=sum(pararrayfun(par.parallel,@calcP1a,para),2);
		else
			ptotn=sum(pararrayfun(par.parallel,@calcP1ao,para),2);
		endif;
	endif;
endfunction;

function press=calcP1o(tpd,parv)
	[T,theta,p0]=modelTPD1(tpd.T,parv(2),parv(1),parv(4),tpd.rate);
	press=p0*parv(3);
endfunction;

function press=calcP1ao(para)
	[T,theta,p0]=modelTPD1(para.T,para.b,para.a,para.d,para.rate);
	press=p0*para.c;
endfunction;

function press=calcP1(tpd,parv)
	[T,theta,p0]=modelTPDlsode(tpd.T,parv(2),parv(1),parv(4),tpd.rate);
	press=p0*parv(3);
endfunction;

function press=calcP1a(para)
	[T,theta,p0]=modelTPDlsode(para.T,para.b,para.a,para.d,para.rate);
	press=p0*para.c;
endfunction;

