function dxdt=hhode(t,x)
%HH egyenlet megoldása:

%bemenetek:
n=x(1); % K csatorna aktiválódás ([0,1] arányszám)
m=x(2); % Na csatorna aktiválódás ([0,1] arányszám)
h=x(3); % Na csatorna inaktiválódás ([0,1] arányszám)
V=x(4); % membrán potenciál (mV)

%fix paraméterek (állandók kísérletekbõl):
C_m=1; % membrán kapacitás
% csatornákhoz kapcsolódó maximális vezetõképességek:
gK=36; 
gNa=120;
gCl=0.3;
% egyensúlyi reverz potenciálok az ionokra:
V_K=-82;
V_Na=45;
V_Cl=-59;

dxdt=zeros(4,1);

%kimenetek a HH egyenletek alapján: 
% alphan(V), betan(V), stb... függvények meghívandók az egyenleteken belül (kísérleti paraméterek mögöttük) 
dxdt(1) = ; % dn : K csat aktivitás változás 
dxdt(2) = ; % dm : Na csat aktivitás változás
dxdt(3) = ; % dh : Na csat inaktivitás változás
dxdt(4) = ; %dV : membránpotenciál változás