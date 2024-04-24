function [c,f,s] = pde_diffusion_1(x,t,u,ux) 
%  PDE_DIFFUSION_1 
%% 
% Diffúziós együttható
D = 1;
%% 
% A három függvény értéke, amelyek meghatározzák a PDE alakját:
c = 1;
f = D*ux;
s = 0;