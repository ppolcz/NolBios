% 12. gyakorlat:
% Fitzhugh-Nagumo neuron mmodell (egyszerûsített HH)
% dx=c*(x-x^3/3+y) gyors válasz (HH m és V)
% dy=(1/c)*(a-x-b*y) lassú válasz (HH n és h)
% paraméterek:
a=;
b=;
c=;

% Feladatok:
% Implementáld a FN egyenletrendszert
% majd vizsgáld meg az alábbiak szerint: 

%% 
% 1. Indítsd random helyekrõl trajektóriákat az egyensúlyi pont körül:
% Mit jelez a vektormezõ? Hogy látszik rajta a gyors és a lassú komponens?
% Mikor alakul ki akciós potenciál?

% 2. Egyensúly körüli helyekrõl indítok trajektóriákat: kerülje meg a hurkot, mikor van AP? 

% 3. görbét tologatom - Áram komponens hozzáadásával, bifurkációk:
% Adjunk áramot (I) a rendszerbe: dv=c*(x-x^3/3+y+I)
% Nézzük meg hogyan változik a rendszer viselkedése! Hol vannak bifurkációs pontok? 
% Itt az áram változzon a ciklusban.

figure;
 for i=1:20
     
    % dx=0 dy=0 egyeneseket megad:
    vx=-2.5:0.05:2.5;
    vy=; %dx=0 egyenes egyenlete
    wx=vx;
    wy=; %dy=0 egyenes egyenlete

    %egyeneseket ábrázol:
    subplot(2,1,1);
    plot(vx,vy, 'b');
    hold on;
    plot(wx,wy, 'r');
    xlim([-2.5 2.5]);
    ylim([-1.5 1.5]);
    
    %rácspontok megadása:
    xv=-2.5:0.2:2.5;
    xw=xv;
    [] = meshgrid(); %meshgrid
    %a vektormezõ gradiense az adott pontokban:
    
    quiver(); %rácspontokban a vektormezõ ábrázolása

     x0=; %az egyensúlyi pont körüli kezdõpontok
     t0=[0 150];
     [t, phi]=ode45(); %megoldani az adott kezdõpontból az egyenletet
     
     % (x,y) trajektória ábrázolása 
     plot(phi(:,1),phi(:,2),'r'); %ábrázol megoldás trajektóriát
     plot(x0(1), x0(2), 'mx') %kezdõpont
     title(['kezdõpont: ', num2str(x0)]);
     xlabel('gyors válasz (v)'),
     ylabel('lassú válasz (w)');
     
     subplot(2,1,2);
     % a gyors válasz (~V változása, maga az AP) ábrázolása
     hold on
     plot(t,phi(:,1)); % phi(:,1) a v-hez tartozó gyors változás az idõ fvnyében
     ylim([-4 2.5]);
     xlabel('idõ'),
     ylabel('gyors válasz (v)');
     
     waitforbuttonpress;
 end

%%
% 4. Hiperpolarizáció: dv=c*(x-x^3/3+y-I) negatív áram alatt tartom a rendszert egy ideig 
% Készítsük el az FNode függvénybe a FN egyenletnek azt a változatát, ahol
% adott "hip" ideig "Ihip" mértékû árammal "lefojtjuk" a sejtet (ez a
% hiperpolarizáció), majd ez után az idõ után megszüntetjük a negatív
% hatást. Mi történik ekkor?

figure;
    hip=; %ennyi ideig fojtom le a sejtet
    Ihip=; %ekkora árammal
    t0=[0 150];
    x0=[-2 0.625-1];
    [t, phi]=ode45(@(t,x)FNode(t,x, hip,Ihip, a,b,c),t0,x0); %ode45 változói (t,x), konstansok (hip, Ihip, a,b,c) 
    
    for ti=1:5:length(t) %adott ideig ábrázolom a trajektóriát, mindig kicsit tovább, így lesz folyamatos az animáció
        if t(ti)>hip %normál eset, normál egynsúlyi pont
            vx=-2.5:0.05:2.5;
            vy=;
            wx=vx;
            wy=;
            subplot(2,1,1);
            plot(vx,vy,'b',wx,wy,'b');            
        else %hiperpolarizált est, hiperpolarizált egyensúlyi pont
             vx=-2.5:0.05:2.5;
             vy=;
             wx=vx;
             wy=;
             subplot(2,1,1);
             plot(vx,vy,'g', wx, wy, 'g');
        end
        hold on;
        plot(phi(1:ti,1),phi(1:ti,2),'r');
        title(['kezdõpont: ', num2str(x0)]);
        xlabel('gyors válasz (v)'),
        ylabel('lassú válasz (w)');   
        
        subplot(2,1,2);
        hold on
        plot(t(1:ti),phi(1:ti,1),'r');
        xlim([0 150]);
        ylim([-4 2.5]);
        xlabel('idõ');
        ylabel('gyors válasz (v)'),
         
        pause(0.05)        
    end


