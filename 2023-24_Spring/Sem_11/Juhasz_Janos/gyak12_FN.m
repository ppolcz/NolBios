% 12. gyakorlat:
% Fitzhugh-Nagumo neuron mmodell (egyszer�s�tett HH)
% dx=c*(x-x^3/3+y) gyors v�lasz (HH m �s V)
% dy=(1/c)*(a-x-b*y) lass� v�lasz (HH n �s h)
% param�terek:
a=-0.7; %-1,0,1
b=0.8;
c=12.5;

% Feladatok:
% Implement�ld a FN egyenletrendszert
% majd vizsg�ld meg az al�bbiak szerint: 

%% 
% 1. Ind�tsd random helyekr�l trajekt�ri�kat az egyens�ly k�r�l:
% Mit jelez a vektormez�? Hogy l�tszik rajta a gyors �s a lass� komponens?
% A gyors komponens ment�n sokkal er�sebb a vektormez� (szinte v�zszintes gradiensek)
% Mikor alakul ki akci�s potenci�l?
% 1. benyom�s:
% Ha a harmadrend� g�rbe f�l�tt van a kezd�pont, akkor "fordul �t" a p�lya,
% -> akkor van kit�r�s a gyors komponensben (ami a fesz�lts�g v�ltoz�st is tartalmazza)
% -> teh�t ekkor van AP 

figure;
for i = 1:20
     
    % dx=0 dy=0 egyeneseket megad:
    vx=-2.5:0.05:2.5;
    vy=-vx+vx.^3/3; %dx=0 egyenes egyenlete
    wx=vx;
    wy=(-wx+a)./b; %dy=0 egyenes egyenlete

    %egyeneseket �br�zol:
    subplot(2,1,1);
    plot(vx,vy, 'b');
    hold on;
    plot(wx,wy, 'b');
    xlim([-2.5 2.5]);
    ylim([-1.5 1.5]);
    
    %r�cspontok megad�sa:
    xv=-2.5:0.2:2.5;
    xw=xv;
    [v1, w2] = meshgrid(xv, xw); %meshgrid
    %a vektormez� gradiense az adott pontokban:
    v1dot=c*(v1-v1.^3/3+w2);
    w2dot=(1/c)*(a-v1-b*w2);
    quiver(v1,w2,v1dot, w2dot,2, 'g'); %r�cspontokban a vektormez� �br�zol�sa

     x0=[-1.2 0.625]+0.2*randn(1,2); %az egyens�lyi pont k�r�li kezd�pontok
     t0=[0 150];
     [t, phi]=ode45(@(t,x) [c*(x(1)-x(1)^3/3+x(2)); (1/c)*(a-x(1)-b*x(2))],t0,x0); %megoldani az adott kezd�pontb�l az egyenletet
     
     % (x,y) trajekt�ria �br�zol�sa 
     plot(phi(:,1),phi(:,2),'r'); %�br�zol megold�s trajekt�ri�t
     plot(x0(1), x0(2), 'mx') %kezd�pont
     title(['kezd�pont: ', num2str(x0)]);
     xlabel('gyors v�lasz (v)'),
     ylabel('lass� v�lasz (w)');
    
     subplot(2,1,2);
     % a gyors v�lasz (~V v�ltoz�sa, maga az AP) �br�zol�sa
     hold on
     plot(t,phi(:,1)); % phi(:,1) a v-hez tartoz� gyors v�ltoz�s az id� fvny�ben
     %plot(t,phi(:,2)); % phi(:,2) a w-hez tartoz� lass� v�ltoz�s az id� fvny�ben
     ylim([-4 2.5]);
     xlabel('id�'),
     ylabel('gyors v�lasz (v)');
     %ylabel('lass� v�lasz (w)');
     
     waitforbuttonpress;
    
 end

%%
% 2. Egyens�ly k�r�li helyekr�l ind�tok trajekt�ri�kat: ker�lje meg a hurkot, mikor van AP? 
% az 1. tipp�nkn�l komplexebb a k�p, a 3rend� g�rbe k�rny�k�n m�s a lass� v�ltoz�nak is van hat�sa...  

figure;  
 for i=1:20
     x0=[-1.2 0.625+0.01*i];%[-1.2 0.6903+0.00001*i];%
     t0=[0 50];
      [t, phi]=ode45(@(t,x) [c*(x(1)-x(1)^3/3+x(2)); (1/c)*(a-x(1)-b*x(2))],t0,x0); %megold kezd�pontb�l egyenletet
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
     xlabel('gyors v�lasz (v)'),
     ylabel('lass� v�lasz (w)');
     title(['kezd�pont: ', num2str(x0)]);
     
     
     subplot(2,1,2);
     xlim([0 50]);
     ylim([-4 2.5]);
     hold on
     plot(t,phi(:,1));
     xlabel('id�'),
     ylabel('gyors v�lasz (v)');
     
     waitforbuttonpress;
 end
 
%%
% 3. g�rb�t tologatom - �ram komponens hozz�ad�s�val, bifurk�ci�k:
% Adjunk �ramot (I) a rendszerbe: dv=c*(x-x^3/3+y+I)
% N�zz�k meg hogyan v�ltozik a rendszer viselked�se! Hol vannak bifurk�ci�s pontok? 
% stabil pont -> Hopf bifurk�ci� (p�lya sz�letik) -> stabil, vonz� peridikus p�lya -> Hopf bifurk�ci� (p�lya meghal) -> stabil pont 
% Hopf bifurk�ci�: �tmenet egyens�lyi pont �s periodikus megold�s k�zt

figure;
for I=0:0.1:2 % �ram mennyis�ge
 
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
 title(['�ram (eltol�s): ', num2str(I)]);
 
 x0=[-1 -1];
 t0=[0 150];
 
 [t, phi]=ode45(@(t,x) [
        c*(x(1)-x(1)^3/3+x(2)+I)
        (1/c)*(a-x(1)-b*x(2))
        ],t0,x0);
     subplot(2,1,1);
     %hold on
     plot(phi(:,1),phi(:,2),'r');
     xlabel('gyors v�lasz (v)'),
     ylabel('lass� v�lasz (w)');
     hold off;
     
     subplot(2,1,2);
     %hold on
     plot(t,phi(:,1));
     ylim([-4 4]);
     xlabel('id�');
     ylabel('gyors v�lasz (v)'),
     
     waitforbuttonpress;
end

%%
% 4. Hiperpolariz�ci�: dv=c*(x-x^3/3+y-I) negat�v �ram alatt tartom a rendszert egy ideig 
% K�sz�ts�k el az FNode f�ggv�nybe a FN egyenletnek azt a v�ltozat�t, ahol
% adott "hip" ideig "Ihip" m�rt�k� �rammal "lefojtjuk" a sejtet (ez a
% hiperpolariz�ci�), majd ez ut�n az id� ut�n megsz�ntetj�k a negat�v
% hat�st. Mi t�rt�nik ekkor?
% hiperpolariz�st sejt nem t�zel, ha ennek v�ge, akkor m�r l�trej�het az AP
% (megfelel� param�terek eset�n) k�ls� gerjeszt�s n�lk�l is
% az egyens�lyi pontokat v�ltoztattuk
% gyors, lass� ir�ny is l�tszik
figure;
    x1=-2;
    hip=50; %ennyi ideig fojtom le a sejtet
    Ihip=-0.6; %ekkora �rammal
    t0=[0 150];
    x0=[x1 0.625-1];
    [t, phi]=ode45(@(t,x)FNode_kesz(t,x, hip,Ihip, a,b,c),t0,x0); %ode45 v�ltoz�i (t,x), konstansok (hip, Ihip, a,b,c) 
    
    for ti=1:5:length(t) %adott ideig �br�zolom a trajekt�ri�t, mindig kicsit tov�bb, �gy lesz folyamatos az anim�ci�
        if t(ti)>hip %norm�l eset, norm�l egyns�lyi pont
            vx=-2.5:0.05:2.5;
            vy=-vx+vx.^3/3;
            wx=vx;
            wy=(a-wx)./b;
            subplot(2,1,1);
            plot(vx,vy,'b',wx,wy,'b');            
        else %hiperpolariz�lt est, hiperpolariz�lt egyens�lyi pont
             vx=-2.5:0.05:2.5;
             vy=-vx+vx.^3/3 - Ihip ;
             wx=vx;
             wy=(a-wx)./b;
             subplot(2,1,1);
             plot(vx,vy,'g', wx, wy, 'g');
        end
        hold on;
        plot(phi(1:ti,1),phi(1:ti,2),'r');
        title(['kezd�pont: ', num2str(x0)]);
        xlabel('gyors v�lasz (v)'),
        ylabel('lass� v�lasz (w)');   
        
        subplot(2,1,2);
        hold on
        plot(t(1:ti),phi(1:ti,1),'r');
        xlim([0 150]);
        ylim([-4 2.5]);
        xlabel('id�');
        ylabel('gyors v�lasz (v)'),
         
        pause(0.05)        
    end


