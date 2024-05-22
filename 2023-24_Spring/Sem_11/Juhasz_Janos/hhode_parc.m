function dxdt=hhode_parc(t,x)
%HH egyenlet megol�sa:

d = 1; %t�rbeli l�p�sk�z
%fix param�terek (�lland�k k�s�rletekb�l):
C_m=1; % membr�n kapacit�s
% csatorn�khoz kapcsol�d� maxim�lis vezet�k�pess�gek:
gk=36;
gNa=120;
gCl=0.3;
% egyens�lyi reverz potenci�lok az ionokra:
V_k=-82;
V_Na=45;
V_Cl=-59;
ujinger_t=20; % ekkor j�n a 2. impulzus
ujinger=-20; % ekkora az �j gerjeszt� inger (potenci�l (mV)) az axon elej�n
%ujinger_t=10 id�pillanat m�lva; V=-20 -> nincs AP 
%ujinger_t=20 ; V=-20 -> van AP 
%ujinger_t=10 ; V=-10 -> van AP (nagyobb �rammal ingerelhet� m�r)
%ujinger_t=5 ; V=-10 -> nics AP, hiperpolariz�lt m�g a membr�n 

% inger�let terjed�se v�gis az axonon: HH egyenletek, a dV egyenletben plusz diff�zi�s tag a potenci�l terjed�s��rt 
% (a D diff�zi�s �lland�t �s a dt-t most vegy�k 1-nek) 

%
%
%adott id�ben (ujinger_t), adott m�ret� ingerl�s (V) az axon elj�re (x(4))
if abs(t-ujinger_t)<0.1
x(4)=ujinger;
end
%}
%kicsit bug-os n�hol...
x=reshape(x, 4, []); %bemenetek m�trix� m�retez�se (minden n az 1. sorban, minden m a 2. sorban, stb)
dxdt=zeros(size(x));
n=x(1,:); % K csatorna aktiv�l�d�s minden pontban ([0,1] ar�nysz�m)
m=x(2,:); % Na csatorna aktiv�l�d�s minden pontban ([0,1] ar�nysz�m)
h=x(3,:); % Na csatorna inaktiv�l�d�s minden pontban ([0,1] ar�nysz�m)
V=x(4,:); % membr�n potenci�l (mV)
V_before=[V(1),V(1:end-1)]; % zero flux hat�r
V_after=[V(2:end),V(end)]; % zero flux hat�r
%kimenetek a HH egyenletek alapj�n: 
% alphan(V), betan(V), stb... f�ggv�nyek megh�vand�k az egyenleteken bel�l (k�s�rleti param�terek m�g�tt�k) 
dxdt(1,:) = ; % dn : K csat aktivit�s v�ltoz�s minden pontban 
dxdt(2,:) = ; % dm : Na csat aktivit�s v�ltoz�s minden pontban
dxdt(3,:) = ; % dh : Na csat inaktivit�s v�ltoz�s minden pontban
dxdt(4,:) = ; %dV : membr�npotenci�l v�ltoz�s + 1D-s diff�zi�
dxdt=reshape(dxdt,[],1); %kimenet oszlopp� visszaalak�t�sa
%}

%{
% cs�ny�bb, de kev�sb� bug-os

dxdt=zeros(length(x),1);

% axon eleje:
n=x(1);
m=x(2);
h=x(3);
V=x(4);

%{
%adott id�ben (ujinger_t), adott m�ret� ingerl�s (V) az axon elj�re (x(4))
if abs(t-ujinger_t)<0.1
V=ujinger;
end
%}

dxdt(1) = ;
dxdt(2) = ;
dxdt(3) = ;
dxdt(4) = ;

% axon k�ztes r�sze
for i = 1:length(x)/4-2
    n= ;
    m= ;
    h= ;
    V= ;
    dxdt(4*i+1) = ;
    dxdt(4*i+2) = ;
    dxdt(4*i+3) = ;
    dxdt(4*i+4) = ;
    
end

% axon v�ge
n=x(end-3);
m=x(end-2);
h=x(end-1);
V=x(end);
dxdt(end-3) = ;
dxdt(end-2) = ;
dxdt(end-1) = ;
dxdt(end) = ;
%}
