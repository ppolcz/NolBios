% 12. gyakorlat:
% Fitzhugh-Nagumo neuron mmodell (egyszerûsített HH)
% dx=c*(x-x^3/3+y) gyors válasz (HH m és V)
% dy=(1/c)*(a-x-b*y) lassú válasz (HH n és h)
% paraméterek:
a=-0.7; %-1,0,1
b=0.8;
c=12.5;

% Feladatok:
% Implementáld a FN egyenletrendszert
% majd vizsgáld meg az alábbiak szerint: 

%% 
% 1. Indítsd random helyekrõl trajektóriákat az egyensúly körül:
% Mit jelez a vektormezõ? Hogy látszik rajta a gyors és a lassú komponens?
% A gyors komponens mentén sokkal erõsebb a vektormezõ (szinte vízszintes gradiensek)
% Mikor alakul ki akciós potenciál?
% 1. benyomás:
% Ha a harmadrendû görbe fölött van a kezdõpont, akkor "fordul át" a pálya,
% -> akkor van kitérés a gyors komponensben (ami a feszültség változást is tartalmazza)
% -> tehát ekkor van AP 

figure;
for i = 1:20
     
    % dx=0 dy=0 egyeneseket megad:
    vx=-2.5:0.05:2.5;
    vy=-vx+vx.^3/3; %dx=0 egyenes egyenlete
    wx=vx;
    wy=(-wx+a)./b; %dy=0 egyenes egyenlete

    %egyeneseket ábrázol:
    subplot(2,1,1);
    plot(vx,vy, 'b');
    hold on;
    plot(wx,wy, 'b');
    xlim([-2.5 2.5]);
    ylim([-1.5 1.5]);
    
    %rácspontok megadása:
    xv=-2.5:0.2:2.5;
    xw=xv;
    [v1, w2] = meshgrid(xv, xw); %meshgrid
    %a vektormezõ gradiense az adott pontokban:
    v1dot=c*(v1-v1.^3/3+w2);
    w2dot=(1/c)*(a-v1-b*w2);
    quiver(v1,w2,v1dot, w2dot,2, 'g'); %rácspontokban a vektormezõ ábrázolása

     x0=[-1.2 0.625]+0.2*randn(1,2); %az egyensúlyi pont körüli kezdõpontok
     t0=[0 150];
     [t, phi]=ode45(@(t,x) [c*(x(1)-x(1)^3/3+x(2)); (1/c)*(a-x(1)-b*x(2))],t0,x0); %megoldani az adott kezdõpontból az egyenletet
     
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
     %plot(t,phi(:,2)); % phi(:,2) a w-hez tartozó lassú változás az idõ fvnyében
     ylim([-4 2.5]);
     xlabel('idõ'),
     ylabel('gyors válasz (v)');
     %ylabel('lassú válasz (w)');
     
     waitforbuttonpress;
    
 end

%%
% 2. Egyensúly körüli helyekrõl indítok trajektóriákat: kerülje meg a hurkot, mikor van AP? 
% az 1. tippünknél komplexebb a kép, a 3rendû görbe környékén más a lassú változónak is van hatása...  

figure;  
 for i=1:20
     x0=[-1.2 0.625+0.01*i];%[-1.2 0.6903+0.00001*i];%
     t0=[0 50];
      [t, phi]=ode45(@(t,x) [c*(x(1)-x(1)^3/3+x(2)); (1/c)*(a-x(1)-b*x(2))],t0,x0); %megold kezdõpontból egyenletet
     %xlim([1.185 1.202]);
     %ylim([-1 1]);
     
     subplot(2,1,1);
     vx=-2.5:0.05:2.5;
     vy=-vx+vx.^3/3;
     wx=vx;
     wy=(a-wx)./b;
     plot(vx,vy,'-b', wx,wy, '-b');
     
     hold on;
     plot(phi(:,1),phi(:,2),'r', x0(1), x0(2), 'mx');
     xlabel('gyors válasz (v)'),
     ylabel('lassú válasz (w)');
     title(['kezdõpont: ', num2str(x0)]);
     
     
     subplot(2,1,2);
     xlim([0 50]);
     ylim([-4 2.5]);
     hold on
     plot(t,phi(:,1));
     xlabel('idõ'),
     ylabel('gyors válasz (v)');
     
     waitforbuttonpress;
 end
 
%%
% 3. görbét tologatom - Áram komponens hozzáadásával, bifurkációk:
% Adjunk áramot (I) a rendszerbe: dv=c*(x-x^3/3+y+I)
% Nézzük meg hogyan változik a rendszer viselkedése! Hol vannak bifurkációs pontok? 
% stabil pont -> Hopf bifurkáció (pálya születik) -> stabil, vonzó peridikus pálya -> Hopf bifurkáció (pálya meghal) -> stabil pont 
% Hopf bifurkáció: átmenet egyensúlyi pont és periodikus megoldás közt

figure;
for I=0:0.1:2 % áram mennyisége
 
 vx=-2.5:0.05:2.5;
 vy=-vx+vx.^3/3-I;
 wx=vx;
 wy=(a-wx)./b;
 
 subplot(2,1,1);
 plot(vx,vy);
 hold on;
 plot(wx,wy);
 xlim([-2.5 2.5]);
 ylim([-6 6]);
 title(['Áram (eltolás): ', num2str(I)]);
 
 x0=[-1 -1];
 t0=[0 150];
 
 [t, phi]=ode45(@(t,x) [
        c*(x(1)-x(1)^3/3+x(2)+I)
        (1/c)*(a-x(1)-b*x(2))
        ],t0,x0);
     subplot(2,1,1);
     %hold on
     plot(phi(:,1),phi(:,2),'r');
     xlabel('gyors válasz (v)'),
     ylabel('lassú válasz (w)');
     hold off;
     
     subplot(2,1,2);
     %hold on
     plot(t,phi(:,1));
     ylim([-4 4]);
     xlabel('idõ');
     ylabel('gyors válasz (v)'),
     
     waitforbuttonpress;
end

%%
% 4. Hiperpolarizáció: dv=c*(x-x^3/3+y-I) negatív áram alatt tartom a rendszert egy ideig 
% Készítsük el az FNode függvénybe a FN egyenletnek azt a változatát, ahol
% adott "hip" ideig "Ihip" mértékû árammal "lefojtjuk" a sejtet (ez a
% hiperpolarizáció), majd ez után az idõ után megszüntetjük a negatív
% hatást. Mi történik ekkor?
% hiperpolarizást sejt nem tüzel, ha ennek vége, akkor már létrejöhet az AP
% (megfelelõ paraméterek esetén) külsõ gerjesztés nélkül is
% az egyensúlyi pontokat változtattuk
% gyors, lassú irány is látszik
figure;
    x1=-2;
    hip=50; %ennyi ideig fojtom le a sejtet
    Ihip=-0.6; %ekkora árammal
    t0=[0 150];
    x0=[x1 0.625-1];
    [t, phi]=ode45(@(t,x)FNode_kesz(t,x, hip,Ihip, a,b,c),t0,x0); %ode45 változói (t,x), konstansok (hip, Ihip, a,b,c) 
    
    for ti=1:5:length(t) %adott ideig ábrázolom a trajektóriát, mindig kicsit tovább, így lesz folyamatos az animáció
        if t(ti)>hip %normál eset, normál egynsúlyi pont
            vx=-2.5:0.05:2.5;
            vy=-vx+vx.^3/3;
            wx=vx;
            wy=(a-wx)./b;
            subplot(2,1,1);
            plot(vx,vy,'b',wx,wy,'b');            
        else %hiperpolarizált est, hiperpolarizált egyensúlyi pont
             vx=-2.5:0.05:2.5;
             vy=-vx+vx.^3/3 - Ihip ;
             wx=vx;
             wy=(a-wx)./b;
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


