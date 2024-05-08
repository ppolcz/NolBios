function [G] = RegularNetwork(rajz, fokszam, nv, k, p)
%Példa fv hívások: [G] = RegularNetwork(1, 0, 10, 2, 0.2);
%[G] = RegularNetwork(0, 1, 1500, 2, 0.2) %fokszám eloszlás Poisson, egyre szélesebb minnél nagyobb az áthuzalozás

%   RegularNetwork metódus gráfgenerálásra használható
%
%   Bemeneti paraméterek:
%   nv      - csúcsok száma (vertices)
%   k       - szomszédság
%   p       - átvezetékezés valószínûsége
%   rajz    - 1: kirajzolás lesz; 0: nem
%   fokszam - 1: fokszámok eloszlását számol és rajzol; 0: nem
%
%   Kimeneti paraméter egy G adatstruktúra a következõ elemekkel:
%   G.Adj   - a gráf szomszédsági mátrixa (1 - a két csúcs között van él, 0
%   - a két csúcs között nincs él
%   G.nv    - csúcsok száma (vertices)
%   G.ne    - élek száma (edges)

% létrehozzuk az alapmátrixot csupa nullával
A=zeros(nv,nv);

% k, a fokszám hibakezelése
k=max(abs(k),2);
if mod(k,2)==1
    k=k+1;
end

for k_=1:k
    % hozzáadjuk az adott csúcshoz tartozó kis diagonálmátrixot
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

%p=0.2; %átvezetékezés valószínûsége
if nargin==5
[G] = WSmodel(rajz, fokszam, G, p);
end
