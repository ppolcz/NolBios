%9. gyakorlat
% SI rendszer egyenlet utazóhulláma:
dt = 0.1; %idõbeli lépésköz 
dx = 1; %térbeli lépésköz

x = 1:dx:1000; %viszgált pontok az egyenesen
nx = length(x);
tmax = 150;
nt = round(tmax/dt); %idõtartomány vége
indv = zeros(nt,1); %sebességhez adatok ide;

c = 1; %hulámterjedés sebességei állandó
D = 4; %diffúzió sebessége
tau = 2; %fertõzési ráta

I = zeros(nx,1); %a betegek kezdeti értékei
%kezdeti értékek:
I(1) = 0.1; %kezdeti fertõzöttek a tér elején
I(250:300) = 0.001; %1;%%0.5; %0.1; %0.01;% % %kis területrõl indul %fertõzöttek kezdeti értéktõl nem függ a hullámfront 
%eta(1:100) = linspace(0.3,0,100); %a kezdeti hulámfronttól se függ a hullám alak
S = 1-I; %fertõzhetõk kezdeti értékei
%fi_kezd = fi(1);

figure;
%kezdeti érték rajzol:
plot(1:nx, I, 'r', 1:nx, S, 'g', 'LineWidth', 2);
grid(gca, 'minor');
    ylim([-0.2 1.2]);
    xlim([1 1000]);
    title('piros: fertõzött (\eta), zöld: fertõzhetõ (\phi); t=0');
    xlabel('x: 1D tér');
    waitforbuttonpress;
    
for k = 1:nt %minden idõpillanatban
    S_kp1 = S;
    I_kp1 = I;
    for i = 1:nx %minden pontban kiszámolom a hullám aktuális értékét (a széleket külön kell kezelni)
        if i==1
            elott=1;
            utan=i+1;
        elseif i==nx 
            elott=i-1;
            utan=nx;
        else
            elott=i-1;
            utan=i+1;
        end
        % dfi/dt=c*fi*(-tau*eta)+D*dfi^2/dx^2 %fertõzhetõek (S) változása, reakció(-betegedés)+diffúzió  
        % deta/dt=c*eta*(tau*fi-1)+D*deta^2/dx^2 %fertõzöttek (I) változása, reakció(+betegedés-gyógyulás)+diffúzió 
        S_kp1(i) = S(i) +c*dt*S(i)*(-tau*I(i))+D*dt/dx^2*(S(utan)-2*S(i)+S(elott));
        I_kp1(i) = I(i) +c*dt*I(i)*(tau*S(i)-1)+D*dt/dx^2*(I(utan)-2*I(i)+I(elott));
    end
    S = S_kp1;
    I = I_kp1;
    
    %ábrázol:
    plot(1:nx, I, 'r', 1:nx, S, 'g', 'LineWidth', 2);
    grid(gca, 'minor');
    ylim([-0.2 1.2]);
    xlim([1 1000]);
    title(['piros: fertõzött (\eta), zöld: fertõzhetõ (\phi); t=', num2str(k*dt),]);
    xlabel('x: 1D tér');
    %waitforbuttonpress;
     pause(0.01);
     
     %{
    % hullám sebességét meghatároz: 
 kul=abs(fi-0.5); 
 [mi, hol]=min(kul); %megkeresem a hullám közepét (ami legközelebb van 0.5-höz)
 indv(t,1)=hol(1); %tárolom, hogy hol, milyen x-nél van ez adott idõpillanatban
    %}
end
 
%{
%sebesség:
elvart_min_sebesseg=2*sqrt(tau*fi_kezd-1) %u(1)=1; többi 0 kezdeti esetben
 tkezd=1/4*nt; %egy intervallumon belül, ahol stabil a hullám sebessége
 tveg=4/5*nt;
sebesseg=abs(indv(tveg,1)-indv(tkezd,1))/abs((tveg-tkezd)*dt) %mennyit halad a hullám közepe adott idõ alatt 
%}




