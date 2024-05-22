function dxdt=hhode_parc(t,x)
%HH egyenlet megolása:

d = 1; %térbeli lépésköz
%fix paraméterek (állandók kísérletekbõl):
C_m=1; % membrán kapacitás
% csatornákhoz kapcsolódó maximális vezetõképességek:
gk=36;
gNa=120;
gCl=0.3;
% egyensúlyi reverz potenciálok az ionokra:
V_k=-82;
V_Na=45;
V_Cl=-59;
ujinger_t=20; % ekkor jön a 2. impulzus
ujinger=-20; % ekkora az új gerjesztõ inger (potenciál (mV)) az axon elején
%ujinger_t=10 idõpillanat múlva; V=-20 -> nincs AP 
%ujinger_t=20 ; V=-20 -> van AP 
%ujinger_t=10 ; V=-10 -> van AP (nagyobb árammal ingerelhetõ már)
%ujinger_t=5 ; V=-10 -> nics AP, hiperpolarizált még a membrán 

% ingerület terjedése végis az axonon: HH egyenletek, a dV egyenletben plusz diffúziós tag a potenciál terjedéséért 
% (a D diffúziós állandót és a dt-t most vegyük 1-nek) 

%
%
%adott idõben (ujinger_t), adott méretû ingerlés (V) az axon eljére (x(4))
if abs(t-ujinger_t)<0.1
x(4)=ujinger;
end
%}
%kicsit bug-os néhol...
x=reshape(x, 4, []); %bemenetek mátrixá méretezése (minden n az 1. sorban, minden m a 2. sorban, stb)
dxdt=zeros(size(x));
n=x(1,:); % K csatorna aktiválódás minden pontban ([0,1] arányszám)
m=x(2,:); % Na csatorna aktiválódás minden pontban ([0,1] arányszám)
h=x(3,:); % Na csatorna inaktiválódás minden pontban ([0,1] arányszám)
V=x(4,:); % membrán potenciál (mV)
V_before=[V(1),V(1:end-1)]; % zero flux határ
V_after=[V(2:end),V(end)]; % zero flux határ
%kimenetek a HH egyenletek alapján: 
% alphan(V), betan(V), stb... függvények meghívandók az egyenleteken belül (kísérleti paraméterek mögöttük) 
dxdt(1,:) = ; % dn : K csat aktivitás változás minden pontban 
dxdt(2,:) = ; % dm : Na csat aktivitás változás minden pontban
dxdt(3,:) = ; % dh : Na csat inaktivitás változás minden pontban
dxdt(4,:) = ; %dV : membránpotenciál változás + 1D-s diffúzió
dxdt=reshape(dxdt,[],1); %kimenet oszloppá visszaalakítása
%}

%{
% csúnyább, de kevésbé bug-os

dxdt=zeros(length(x),1);

% axon eleje:
n=x(1);
m=x(2);
h=x(3);
V=x(4);

%{
%adott idõben (ujinger_t), adott méretû ingerlés (V) az axon eljére (x(4))
if abs(t-ujinger_t)<0.1
V=ujinger;
end
%}

dxdt(1) = ;
dxdt(2) = ;
dxdt(3) = ;
dxdt(4) = ;

% axon köztes része
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

% axon vége
n=x(end-3);
m=x(end-2);
h=x(end-1);
V=x(end);
dxdt(end-3) = ;
dxdt(end-2) = ;
dxdt(end-1) = ;
dxdt(end) = ;
%}
