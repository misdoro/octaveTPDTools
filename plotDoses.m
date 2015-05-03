#Plot doses from struct.
#Uses inp.doses, inp.integrals arrays
function plotDoses(inp);
inp.doses*=1e10;
inp.integrals*=1e10;
collage=polyfit(inp.doses,inp.integrals,1);
points=[0;max(inp.doses)];
plot(inp.doses,inp.integrals,'o',points,polyval(collage,points));

np=length(inp.doses);
for idx=1:np
	text(inp.doses(idx),inp.integrals(idx),sprintf(" <%d",np+1-idx));
end

printf("Itpd/Idose=%.2f\n",collage(1));
printf("offset=%.2f\n",collage(2));
xlabel("Dose integrals, ×10-10 AS");
ylabel("TPD integrals, ×10-10 AS");
endfunction;