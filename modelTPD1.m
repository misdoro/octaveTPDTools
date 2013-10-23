function [ Tode, theta, p ] = modelTPD1(T,theta0,v,E)
% Model TPD curve
%
% T - temperatures needed, independent variable
% N - normalization factor, area under p(T)
% v - frequency
% E - energy level in electron volts


k = 1.38e-23*6.2e18; % (J K^-1)*(eV/J) = eV K^-1 = 8.6e-5 eV K^-1
odepar.ek1=E/k;

a = 3/60; % K/s, heating rate
odepar.nuoa=v/a;

odepar.myzero=10*eps;


vopt=odeset("NormControl","on","Stats","off");


tic();
ret=ode45(@odemlde,T,theta0,vopt,odepar);
toc()
#ret.stats
Tode=ret.x;
theta=ret.y;
p=-arrayfun(@odemlde,Tode,theta,odepar);

endfunction;

