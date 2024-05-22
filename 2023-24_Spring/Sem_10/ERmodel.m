function [G,A] = ERmodel(nV,nE)
%P�lda fv h�v�sok: [G] = ERmodel(1, 0, 10, 0, 5)
%[G] = ERmodel(0, 1, 1500, 0, 8000) %a foksz�m eloszl�s Poisson, harangg�rbe
%
%   ERmodel met�dus gr�fgener�l�sra haszn�lhat�, Erd�s - R�nyi elj�r�sa
%   szerint
%
%   Bemenete:
%   nV       - cs�csok sz�ma
%   nE       - az �lek sz�ma

idx = sort(randperm(nV*(nV-1)/2,nE));
d = ceil(idx/nV);
i = mod(idx-1,nV)+1;
j = mod(idx+d-1,nV)+1;

A = sparse(i,j,idx*0+1,nV,nV);
A = A + A';

G = graph(A);

