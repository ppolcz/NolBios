function [G] = BAmodel_v2(n,n0,L,L0)
%Példa fv hívások: [G] = BAmodel(1, 0, 50, 3, 1, 3)
%[G] = BAmodel(0, 1, 1500, 3, 1, 3) %a fokszám eloszlás hatványfüggvény
%
%   BAmodel metódus gráfgenerálásra használható, Barabási - Albert eljárása
%   szerint
%
%   Bemenete:
%   n       - csúcsok száma
%   n0      - kezdeti csúcsszám 
%   L       - az új csúcshoz húzott élek száma
%   L0      - az élek száma a gráf létrehozásakor
%   rajz    - 1: kirajzolás lesz; 0: nem
%   fokszam - 1: fokszámok eloszlását számol és rajzol; 0: nem
%
%   Kimenete egy G adatstruktúra a következõ elemekkel:
%   G.Adj   - a gráf szomszédsági mátrixa (1 - a két csúcs között van él, 0
%   - a két csúcs között nincs él
%   G.nv    - csúcsok száma (vertices)
%   G.ne    - élek száma (edges)

%%


%%

A = zeros(n);

% Elején véletlenszerûen feltöltjük a kezdeti gráfot L0 éllel
A0 = zeros(n0);
idx = find(tril(ones(n0),-1));
sigma = randperm(numel(idx),L0);
A0(idx(sigma)) = 1;
A0 = A0 + A0' - diag(A0);
A(1:n0,1:n0) = A0;

% Mindegyik idõpillanatban hozzáadunk egy csúcsot, míg m <= L, ekkor ez a
% csúcs a meglévõ m csúccsal lehet kapcsolatban. Hogy segítsük a preferencia
% szerinti kapcsolódást, az új csúcsot p(A_i,j) = L_j / sum(A_j)
% valószínûséggel köti a már meglévõ csúcshoz.

for t = find(sum(A) == 0) %új csúcs hozzáadása
    
    k=sum(A); % már meglévõ csúcsok fokszámai
    P   = k /sum(k); %normálás
    PP  = cumsum(P); %valószínûségi eloszlás készítése
    j=0; %behúzott élek száma
    k = 0; %rossz helyre tervezett él
    while (j<L && k<L^3)   % k<L^3 csak a végtelen ciklus elkerülése miatt van
        rv = rand(1);
        index = find(PP >= rv, 1); %melyik csúcshoz tartozzon az él
        % ha még nincs él, hozzáadjuk
        if A(t,index)~=1
            A(t,index)=1;
            A(index,t)=1;
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

G = graph(A);


% MyColorMap = [
%     COL.Color_1
%     COL.Color_Red
%     COL.Color_Dark_Green
%     ];
% colormap(gca,MyColorMap)
