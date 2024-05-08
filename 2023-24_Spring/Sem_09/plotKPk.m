function plotKPk(G)
% G gráf k - P(k) függvényének megjelenítése


A = G.Adj;
n = G.nv;

k=sum(A,2); %fokszámok: az egyes csúcsok kapcsolatainak számai

% histogram készítése a valószínûség megjelenítéséhez
x  = (0:n-1)'; %lehetséges fokszámok
Pk = hist(k,x); %mennyi van az adott fokszmámból
Pk = Pk./n; %normálás összes csúcsszámmal

figure;
bar(x, Pk); %fokszám eloszlás kirajzolása, pl.: bar(...)
xlabel('fokszámok (kapcsolatok száma)');
ylabel('csúcsok aránya ekkora fokszámmal');
%xlim([-1,n-1]);

end

