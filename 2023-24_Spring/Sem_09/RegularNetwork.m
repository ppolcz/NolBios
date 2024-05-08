function [G] = RegularNetwork(rajz, fokszam, nv, k, p)
%P�lda fv h�v�sok: [G] = RegularNetwork(1, 0, 10, 2, 0.2);
%[G] = RegularNetwork(0, 1, 1500, 2, 0.2) %foksz�m eloszl�s Poisson, egyre sz�lesebb minn�l nagyobb az �thuzaloz�s

%   RegularNetwork met�dus gr�fgener�l�sra haszn�lhat�
%
%   Bemeneti param�terek:
%   nv      - cs�csok sz�ma (vertices)
%   k       - szomsz�ds�g
%   p       - �tvezet�kez�s val�sz�n�s�ge
%   rajz    - 1: kirajzol�s lesz; 0: nem
%   fokszam - 1: foksz�mok eloszl�s�t sz�mol �s rajzol; 0: nem
%
%   Kimeneti param�ter egy G adatstrukt�ra a k�vetkez� elemekkel:
%   G.Adj   - a gr�f szomsz�ds�gi m�trixa (1 - a k�t cs�cs k�z�tt van �l, 0
%   - a k�t cs�cs k�z�tt nincs �l
%   G.nv    - cs�csok sz�ma (vertices)
%   G.ne    - �lek sz�ma (edges)

% l�trehozzuk az alapm�trixot csupa null�val
A=zeros(nv,nv);

% k, a foksz�m hibakezel�se
k=max(abs(k),2);
if mod(k,2)==1
    k=k+1;
end

for k_=1:k
    % hozz�adjuk az adott cs�cshoz tartoz� kis diagon�lm�trixot
    A = A + diag(ones(1,length(diag(A,k_))), k_) ...
          + diag(ones(1,length(diag(A,nv-k_))), nv-k_);  
end
A = A + A';

G = struct('Adj', A, 'nv', nv, 'ne', nnz(A)/2);

if rajz==1
    plotGraphBasic(G,5,1);
end

if fokszam==1
    plotKPk(G);
end

%p=0.2; %�tvezet�kez�s val�sz�n�s�ge
if nargin==5
[G] = WSmodel(rajz, fokszam, G, p);
end
