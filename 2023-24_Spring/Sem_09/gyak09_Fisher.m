%9. gyakorlat:
% Fisher egyenlet:
dt = 0.1; %idõbeli lépésköz %D*2*dt/dx^2<1 %0.46 körül "szép szõrösödés", nagy ,kis freki ugyanúgy cseng le
dx = 1; %térbeli lépésköz

x = 1:dx:500; %viszgált pontok az egyenesen
nx = length(x);
tmax=150;%1500;%75;%
nt = round(tmax/dt); %idõtartomány vége
indv=zeros(nt,1); %sebességhez adatok ide;

c =1; %hulámterjedés sebességi állandója
D=1; %diffúzió sebessége

u=zeros(nx,1);%exp(-0.1*x);%  %a hullám kezdeti értékei (ha sehol se 0 kezdetben a rendszer, az befolyásolja a hullámfrontot) 
%kezdeti értékek:
u(1) =1; %5;%0.1;%1;% az 1. pontban 1 (a többiben 0) a hullám kezdetben
u(end) = 1;
% u(1:50) =1;%0.1;%5;%1; % (pozitív) kezdettõl függetlenül 1-re áll be a hullám
% u(250:end)=1;% visszafele is terjed
% u(250:260)=5;%5;%0.1;%1; %csúcsból indul 
% u([250:260, 200:230])=1; %0.1;%5;%1;% 2 csúcs indul 
% u(250:300)=rand(1,51)*2; %random kezdet (nem lehet negatív)
% u(250:300)=linspace(0,2,51); %a kezdeti hulámfronttól se függ hullám alakja
% u(1:100)=linspace(2,0,100); %a kezdeti hulámfronttól se függ
% u(1:100)=linspace(0.5,0,100); %a kezdeti hulámfronttól se függ
% négyzetets Fisher:
% u(250:300)=rand(1,51)*2-1; %random kezdet (lehet negatív)
u(250:300)=linspace(-0.5,2,51); %a kezdeti hulámfronttól se függ (negatív rész is)
%u(200:249)=linspace(-1,1,50); u(250:300)=linspace(1,-0.1,51); u(301:350)=linspace(-0.1,1,50);%majdnem -1-re áll be
%u(200:249)=linspace(-1,1,50); u(250:300)=linspace(1,-0.06,51); u(301:350)=linspace(-0.06,1,50);%sokáig úgy tûnik -1-re áll be, aztán ugrik fel 1-re (hosszabb futás): metastabilitás 

figure;
%kezdeti érték kirajzolása:
plot(u); %hullám értékeinek ábrázolása
grid(gca, 'minor');
    ylim([-1.2 5]);
    xlim([1 450]);
    title('t=0');
    xlabel('x: 1D tér');
    waitforbuttonpress;
    
for t = 1:nt %minden idõpillanatban
    u2 = u; %u2=értékek az elõzõ idõplillanatban; u=új értékek
    for i = 1:nx %minden pontban kiszámolom a hullám aktuális értékét (a széleket külön kell kezelni)
        %kezeljük a peremeket zero flux peremként
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
        % 1. dU/dt = c*U*(1-U) + D*dU^2/dx^2 % reakc tag+diff tag=fisher, térben szétterjed (diff),
        % u(i) = u2(i) + c*dt*u2(i)*(1-u2(i)) + D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % 2. dU/dt = c*U*(1-U^2) + D*dU^2/dx^2 % négyzetes fisher (-1-es stabil pont is van)
        u(i) = u2(i) +c*dt*(u2(i))*(1-u2(i)^2)+D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % hasonló rendszerek:
        % 3. dU/dt=D*dU^2/dx^2 %diff tag csak, szétfolyik ,de nem mozdul
        % u(i) = u2(i) + D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % 4. dU/dt=-c*dU/dx %lin convection: térben utazik, nem folyik szét a hullám  és fenn is marad (reakc) 
        % u(i) = u2(i) - c*dt/dx*(u2(i)-u2(elott));
        % 5. dU/dt=-c*dU/dx+D*dU^2/dx^2 %lin convection+diff=burgers, utazik és szétfolyik 
        % u(i) = u2(i) - c*dt/dx*(u2(i)-u2(elott))+D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        
    end
    
    %ábrázol:
    plot(u);
    grid(gca, 'minor');
    ylim([-1.2 5]);
    xlim([1 450]);
    title(['t=' ,num2str(t*dt)]);
    xlabel('x: 1D tér');
    %waitforbuttonpress;
     pause(0.01);
     
     %
 % fishernél hullám sebességét meghatároz: (1-bõl indított hullámnál mûködik jól)
 kul=abs(u-0.5); 
 [mi, hol]=min(kul); %megkeresem a hullám közepét (ami legközelebb van 0.5-höz)
 indv(t,1)=hol(1); %tárolom, hogy hol, milyen x-nél van ez adott idõpillanatban
    %}
end
 
%
%fishernél:
elvart_sebesseg=2*sqrt(c*D) %u(1)=1; többi 0 kezdeti esetben
tkezd=round(1/4*nt); %egy intervallumon belül, ahol stabil a hullám sebessége (,hogy a tranziensek, peremek ne zavarjanak) 
tveg=round(4/5*nt);
sebesseg=(indv(tveg,1)-indv(tkezd,1))/((tveg-tkezd)*dt) %mennyit halad a hullám közepe adott idõ alatt 
%nem lesz pont az elvárt; 2.036-ra 2 helyett kihozható c=1 D=1, dx=1, dt=0.01 esetén
%}




