function result=dispplotqP(mytpd,param,result);
	hold on;
	printf("Use keys a, d to move the pressure, q to exit\n");
	shiftit=-2;
	inp="0";
	while(inp!="q")
		if (inp=="a")
			shiftit--;
		elseif (inp=="d")
			shiftit++;
		endif
		shiftit
		psi=circshift(mytpd.pi,shiftit);
		plot(psi,mytpd.i,"linewidth",2,"color",mytpd.color);
		[maxi,maxidx]=max(psi);
		maxiq=mytpd.i(maxidx);
		text(maxi,maxiq,strcat(num2str(mytpd.idx)));
		inp=kbhit();
	endwhile
endfunction