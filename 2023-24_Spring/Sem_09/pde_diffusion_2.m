function [c,f,s] = pde_diffusion_2(x,t,u,ux) 

% Diffúziós együttható
D = 1;

% A három függvény értéke, amelyek meghatározzák a PDE alakját:
c = 1;
f = ux;
s = u*(1 - u^2) + ux*10;