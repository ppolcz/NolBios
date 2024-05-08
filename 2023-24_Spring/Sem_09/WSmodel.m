function [G] = WSmodel(rajz, fokszam, RN, p)
%   WSmodel metódus gráfgenerálásra használható, Watts-Strogatz eljárása
%   szerint
%
%   Bemenete:
%   RN      - RegularNetwork metódussal generált gráf
%   p       - átvezetékelés valószínûsége
%   rajz    - 1: kirajzolás lesz; 0: nem
%   fokszam - 1: fokszámok eloszlását számol és rajzol; 0: nem
%
%   Kimenete egy G adatstruktúra a következõ elemekkel:
%   G.Adj   - a gráf szomszédsági mátrixa (1 - a két csúcs között van él, 0
%   - a két csúcs között nincs él
%   G.nv    - csúcsok száma (vertices)
%   G.ne    - élek száma (edges)

% átalakítjuk a mátrixot felsõháromszöggé (mivel irányítatlan a gráf, ezért
% ezt nyugodtan megtehetjük)
A = triu(RN.Adj);

nv = RN.nv;

% megkeressük az élekkel összekötött csúcspárokat
[v1,v2] = find(A);

% ezek közül a következõ éleket kell felbontani (p valószínûséggel számolva) 
dis = (rand(length(v1),1)<=p); 
A(v1(dis),v2(dis))=0; % él kitörlése a szomszédsági mátrixból

nDis = sum(dis); % hány darab élet kell visszarakni, mert törölve lett

pairs_to_disconnect=[v1(dis),v2(dis)];

for n=1:nDis
    % kiválasztjuk az egyik csúcsot a kitörölt élbõl
    i=ceil(rand*size(pairs_to_disconnect,1));
    j=randi([1,2],1);
    v_dis_to_rec=pairs_to_disconnect(i,j);
    
    % keresünk egy nem szomszédos csúcsot
    adj=[find(A(:,v_dis_to_rec)) ; find(A(v_dis_to_rec,:))'];
    non_adj=setdiff(1:nv,adj);
    
    % ne legyen hurok
    v_to_rec = v_dis_to_rec;
    while v_to_rec==v_dis_to_rec
        v_to_rec=non_adj(ceil(rand*length(non_adj)));
    end
        
    % végül beírjuk az A mátrixba az új élt
    S = sort([v_dis_to_rec v_to_rec]);
    A(S(1),S(2))=1;
end

A=A+A'; %négyzetes mátrixá visszaalakítás

G = struct('Adj', A, 'nv', nv, 'ne', nnz(A)/2);

if rajz==1
    plotGraphBasic(G,5,1);
end

if fokszam==1
    plotKPk(G);
end
