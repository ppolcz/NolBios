function [G] = ERmodel(rajz, fokszam, n, p, L)
%Példa fv hívások: [G] = ERmodel(1, 0, 10, 0, 5)
%[G] = ERmodel(0, 1, 1500, 0, 8000) %a fokszám eloszlás Poisson, haranggörbe
%
%   ERmodel metódus gráfgenerálásra használható, Erdõs - Rényi eljárása
%   szerint
%
%   Bemenete:
%   n       - csúcsok száma
%   p       - az él létezésének valószínûsége
%   L       - az élek száma (ha megvan adva, akkor az algoritmus ezt veszi
%   figyelembe
%   rajz    - 1: kirajzolás lesz; 0: nem
%   fokszam - 1: fokszámok eloszlását számol és rajzol; 0: nem
%
%   Kimenete egy G adatstruktúra a következõ elemekkel:
%   G.Adj   - a gráf szomszédsági mátrixa (1 - a két csúcs között van él, 0
%   - a két csúcs között nincs él
%   G.nv    - csúcsok száma (vertices)
%   G.ne    - élek száma (edges)


adj=zeros(n); % üres szomszédsági mátrix

switch nargin
    case 3  % ha csak egy paraméter van, az a csúcsok száma
        % 0.5 a p alapértelmezett értéke most
        for i=1:n
            for j=i+1:n
                if rand<=0.5; adj(i,j)=1; adj(j,i)=1; end
            end
        end
        
    case 4 % ha p is defininált, az él p eséllyel jelenik meg
        for i=1:n
            for j=i+1:n
                if rand<=p; adj(i,j)=1; adj(j,i)=1; end
            end
        end
        
    case 5 % a csúcsoknak és az éleknek pontosan megadott számuk lehet.
        while sum(sum(adj))/2 < L
            i=randi(n); j=randi(n);
            if i==j || adj(i,j)>0; continue; end  % ne legyen se dupla él, se hurok.
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
