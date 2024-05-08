function [G] = BAmodel(rajz, fokszam, n, n0, L, L0)
%P�lda fv h�v�sok: [G] = BAmodel(1, 0, 50, 3, 1, 3)
%[G] = BAmodel(0, 1, 1500, 3, 1, 3) %a foksz�m eloszl�s hatv�nyf�ggv�ny
%
%   BAmodel met�dus gr�fgener�l�sra haszn�lhat�, Barab�si - Albert elj�r�sa
%   szerint
%
%   Bemenete:
%   n       - cs�csok sz�ma
%   n0      - kezdeti cs�csz�m 
%   L       - az �j cs�cshoz h�zott �lek sz�ma
%   L0      - az �lek sz�ma a gr�f l�trehoz�sakor
%   rajz    - 1: kirajzol�s lesz; 0: nem
%   fokszam - 1: foksz�mok eloszl�s�t sz�mol �s rajzol; 0: nem
%
%   Kimenete egy G adatstrukt�ra a k�vetkez� elemekkel:
%   G.Adj   - a gr�f szomsz�ds�gi m�trixa (1 - a k�t cs�cs k�z�tt van �l, 0
%   - a k�t cs�cs k�z�tt nincs �l
%   G.nv    - cs�csok sz�ma (vertices)
%   G.ne    - �lek sz�ma (edges)

A = zeros(n);
steps=n-n0; %ennyi cs�cs kell m�g


% Elej�n v�letlenszer�en felt�ltj�k a kezdeti gr�fot L0 �llel
for j=1:L0
    l_ = 1;
    n_ = 1;
    while (A(l_,n_)~=0) || (l_==n_)
        l_ = randi([1,n0],1);
        n_ = randi([1,n0],1);
    end
    A (l_,n_) = 1;
    A (n_,l_) = 1;
end

% Mindegyik id�pillanatban hozz�adunk egy cs�csot, m�g m <= L, ekkor ez a
% cs�cs a megl�v� m cs�ccsal lehet kapcsolatban. Hogy seg�ts�k a preferencia
% szerinti kapcsol�d�st, az �j cs�csot p(A_i,j) = L_j / sum(A_j)
% val�sz�n�s�ggel k�ti a m�r megl�v� cs�cshoz.

for t = 1:steps %�j cs�cs hozz�ad�sa
    
    k=sum(A(1:n0+t-1,:)); %m�r megl�v� cs�csok foksz�mai
    P   = k /sum(k); %norm�l�s
    PP  = cumsum(P); %val�sz�n�s�gi eloszl�s k�sz�t�se
    j=0; %beh�zott �lek sz�ma
    k = 0; %rossz helyre tervezett �l
    while (j<L && k<L^3)   % k<L^3 csak a v�gtelen ciklus elker�l�se miatt van
        rv = rand(1);
        index = find(PP >= rv, 1); %melyik cs�cshoz tartozzon az �l
        % ha m�g nincs �l, hozz�adjuk
        if A(n0+t,index)~=1
            A(n0+t,index)=1;
            A(index,n0+t)=1;
            j=j+1;
        else
            k=k+1; %% error raise
        end
         
    end
    
    if (k>=L^3)
        display('error');
        G = [];
        return;
    end
    
end

G = struct('Adj', A, 'nv', n, 'ne', nnz(A)/2);

if rajz==1
    plotGraphBasic(G,5,1);
end

if fokszam==1
    plotKPk(G);
end
