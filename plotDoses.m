#Plot doses from struct.
#Uses inp.doses, inp.integrals arrays
function plotDoses(inp);

collage=polyfit(inp.doses,inp.integrals,1)
points=[0;max(inp.doses)];
plot(inp.doses,inp.integrals,'o',points,polyval(collage,points));

for idx=1:length(inp.doses)
	text(inp.doses(idx),inp.integrals(idx),strcat(" <",num2str(idx)));
end

xlabel("Dose integrals, A*S");
ylabel("TPD integrals, A*S");
endfunction;