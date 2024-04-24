function [c,f,s] = pde_diffusion_2(x,t,u,ux) 
%% 
% Diffúziós együttható1f1

D = 1;
%% 
% A három függvény értéke, amelyek meghatározzák a PDE alakját:

c = 1;
f = D*ux;
s = u - u^2;