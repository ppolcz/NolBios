function dxdt=hhode(t,x)
%HH egyenlet megold�sa:

%bemenetek:
n=x(1); % K csatorna aktiv�l�d�s ([0,1] ar�nysz�m)
m=x(2); % Na csatorna aktiv�l�d�s ([0,1] ar�nysz�m)
h=x(3); % Na csatorna inaktiv�l�d�s ([0,1] ar�nysz�m)
V=x(4); % membr�n potenci�l (mV)

%fix param�terek (�lland�k k�s�rletekb�l):
C_m=1; % membr�n kapacit�s
% csatorn�khoz kapcsol�d� maxim�lis vezet�k�pess�gek:
gK=36; 
gNa=120;
gCl=0.3;
% egyens�lyi reverz potenci�lok az ionokra:
V_K=-82;
V_Na=45;
V_Cl=-59;

dxdt=zeros(4,1);

%kimenetek a HH egyenletek alapj�n: 
% alphan(V), betan(V), stb... f�ggv�nyek megh�vand�k az egyenleteken bel�l (k�s�rleti param�terek m�g�tt�k) 
dxdt(1) = ; % dn : K csat aktivit�s v�ltoz�s 
dxdt(2) = ; % dm : Na csat aktivit�s v�ltoz�s
dxdt(3) = ; % dh : Na csat inaktivit�s v�ltoz�s
dxdt(4) = ; %dV : membr�npotenci�l v�ltoz�s