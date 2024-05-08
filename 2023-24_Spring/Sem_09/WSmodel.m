function [G] = WSmodel(rajz, fokszam, RN, p)
%   WSmodel met�dus gr�fgener�l�sra haszn�lhat�, Watts-Strogatz elj�r�sa
%   szerint
%
%   Bemenete:
%   RN      - RegularNetwork met�dussal gener�lt gr�f
%   p       - �tvezet�kel�s val�sz�n�s�ge
%   rajz    - 1: kirajzol�s lesz; 0: nem
%   fokszam - 1: foksz�mok eloszl�s�t sz�mol �s rajzol; 0: nem
%
%   Kimenete egy G adatstrukt�ra a k�vetkez� elemekkel:
%   G.Adj   - a gr�f szomsz�ds�gi m�trixa (1 - a k�t cs�cs k�z�tt van �l, 0
%   - a k�t cs�cs k�z�tt nincs �l
%   G.nv    - cs�csok sz�ma (vertices)
%   G.ne    - �lek sz�ma (edges)

% �talak�tjuk a m�trixot fels�h�romsz�gg� (mivel ir�ny�tatlan a gr�f, ez�rt
% ezt nyugodtan megtehetj�k)
A = triu(RN.Adj);

nv = RN.nv;

% megkeress�k az �lekkel �sszek�t�tt cs�csp�rokat
[v1,v2] = find(A);

% ezek k�z�l a k�vetkez� �leket kell felbontani (p val�sz�n�s�ggel sz�molva) 
dis = (rand(length(v1),1)<=p); 
A(v1(dis),v2(dis))=0; % �l kit�rl�se a szomsz�ds�gi m�trixb�l

nDis = sum(dis); % h�ny darab �let kell visszarakni, mert t�r�lve lett

pairs_to_disconnect=[v1(dis),v2(dis)];

for n=1:nDis
    % kiv�lasztjuk az egyik cs�csot a kit�r�lt �lb�l
    i=ceil(rand*size(pairs_to_disconnect,1));
    j=randi([1,2],1);
    v_dis_to_rec=pairs_to_disconnect(i,j);
    
    % keres�nk egy nem szomsz�dos cs�csot
    adj=[find(A(:,v_dis_to_rec)) ; find(A(v_dis_to_rec,:))'];
    non_adj=setdiff(1:nv,adj);
    
    % ne legyen hurok
    v_to_rec = v_dis_to_rec;
    while v_to_rec==v_dis_to_rec
        v_to_rec=non_adj(ceil(rand*length(non_adj)));
    end
        
    % v�g�l be�rjuk az A m�trixba az �j �lt
    S = sort([v_dis_to_rec v_to_rec]);
    A(S(1),S(2))=1;
end

A=A+A'; %n�gyzetes m�trix� visszaalak�t�s

G = struct('Adj', A, 'nv', nv, 'ne', nnz(A)/2);

if rajz==1
    plotGraphBasic(G,5,1);
end

if fokszam==1
    plotKPk(G);
end
