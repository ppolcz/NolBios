function [G] = ERmodel(rajz, fokszam, n, p, L)
%P�lda fv h�v�sok: [G] = ERmodel(1, 0, 10, 0, 5)
%[G] = ERmodel(0, 1, 1500, 0, 8000) %a foksz�m eloszl�s Poisson, harangg�rbe
%
%   ERmodel met�dus gr�fgener�l�sra haszn�lhat�, Erd�s - R�nyi elj�r�sa
%   szerint
%
%   Bemenete:
%   n       - cs�csok sz�ma
%   p       - az �l l�tez�s�nek val�sz�n�s�ge
%   L       - az �lek sz�ma (ha megvan adva, akkor az algoritmus ezt veszi
%   figyelembe
%   rajz    - 1: kirajzol�s lesz; 0: nem
%   fokszam - 1: foksz�mok eloszl�s�t sz�mol �s rajzol; 0: nem
%
%   Kimenete egy G adatstrukt�ra a k�vetkez� elemekkel:
%   G.Adj   - a gr�f szomsz�ds�gi m�trixa (1 - a k�t cs�cs k�z�tt van �l, 0
%   - a k�t cs�cs k�z�tt nincs �l
%   G.nv    - cs�csok sz�ma (vertices)
%   G.ne    - �lek sz�ma (edges)


adj=zeros(n); % �res szomsz�ds�gi m�trix

switch nargin
    case 3  % ha csak egy param�ter van, az a cs�csok sz�ma
        % 0.5 a p alap�rtelmezett �rt�ke most
        for i=1:n
            for j=i+1:n
                if rand<=0.5; adj(i,j)=1; adj(j,i)=1; end
            end
        end
        
    case 4 % ha p is definin�lt, az �l p es�llyel jelenik meg
        for i=1:n
            for j=i+1:n
                if rand<=p; adj(i,j)=1; adj(j,i)=1; end
            end
        end
        
    case 5 % a cs�csoknak �s az �leknek pontosan megadott sz�muk lehet.
        while sum(sum(adj))/2 < L
            i=randi(n); j=randi(n);
            if i==j || adj(i,j)>0; continue; end  % ne legyen se dupla �l, se hurok.
            adj(i,j)=1; adj(j,i)=adj(i,j);
        end
        
end  

G = struct('Adj', adj, 'nv', n, 'ne', nnz(adj)/2);

if rajz==1
    plotGraphBasic(G,5,1);
end

if fokszam==1
    plotKPk(G);
end
