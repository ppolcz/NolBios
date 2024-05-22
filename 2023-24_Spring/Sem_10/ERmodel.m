function [G,A] = ERmodel(nV,nE)
%Példa fv hívások: [G] = ERmodel(1, 0, 10, 0, 5)
%[G] = ERmodel(0, 1, 1500, 0, 8000) %a fokszám eloszlás Poisson, haranggörbe
%
%   ERmodel metódus gráfgenerálásra használható, Erdõs - Rényi eljárása
%   szerint
%
%   Bemenete:
%   nV       - csúcsok száma
%   nE       - az élek száma

idx = sort(randperm(nV*(nV-1)/2,nE));
d = ceil(idx/nV);
i = mod(idx-1,nV)+1;
j = mod(idx+d-1,nV)+1;

A = sparse(i,j,idx*0+1,nV,nV);
A = A + A';

G = graph(A);

