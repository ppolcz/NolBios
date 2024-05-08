function plotKPk(G)
% G gr�f k - P(k) f�ggv�ny�nek megjelen�t�se


A = G.Adj;
n = G.nv;

k=sum(A,2); %foksz�mok: az egyes cs�csok kapcsolatainak sz�mai

% histogram k�sz�t�se a val�sz�n�s�g megjelen�t�s�hez
x  = (0:n-1)'; %lehets�ges foksz�mok
Pk = hist(k,x); %mennyi van az adott fokszm�mb�l
Pk = Pk./n; %norm�l�s �sszes cs�cssz�mmal

figure;
bar(x, Pk); %foksz�m eloszl�s kirajzol�sa, pl.: bar(...)
xlabel('foksz�mok (kapcsolatok sz�ma)');
ylabel('cs�csok ar�nya ekkora foksz�mmal');
%xlim([-1,n-1]);

end

