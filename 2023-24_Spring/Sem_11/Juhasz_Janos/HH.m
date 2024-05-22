% 11. gyakorlat:
% Hodgkin-Huxley modell akci�s potenci�l kialakul�sa (1 pont vizsg�lata az axon membr�non) (ODE rendszer)
% kezdeti (nyugalmi) �rt�kek: ioncsatorna kapuv�ltoz�k (mennyire nyitott ~ ellen�ll�s) :
% [K csatorna nyit�si, Na csatorna nyit�si, Na csatorna inaktiv�l�d�si, kezdeti membr�npotenci�l]

%Feladatok:
% Implement�ld a HH egyenletrendszert a hhode.m f�jlba!
% �ll�tsd be a kezdeti ingert (impulzus) (�s sz�ks�g eset�n a t�bbi benenetet), �gy hogy
% A) akci�s potenci�l (AP) alakuljon ki,
% B) ne legyen AP (marad a nyugalmi potenci�l az input ut�n), hiperpolariz�lt a membr�n
% C) ne legyen AP, mert t�l magas m�r a membr�npotenci�l,
% D) "cs�pd el" az AP kialakul�s�t (azt az ingerl� potenci�lt, ami m�r �pp el�g az AP-hez),
% Mi jellemz� az AP alakj�ra, megjelen�s�re, az ioncsatorn�k aktivit�s�ra?
% Mik a nyugalmi �rt�keik?

figure;
impulzus=; %start impulzus
for i=1:20
ySS=[0.24 0.03 0.76 -65+impulzus]; %nyugalmi potenci�l, innen indul a rendszer
[t,x]=ode45('hhode',[0 50],ySS); %%HH egyenlet (k�l�n f�jl: hhode.m)
subplot(2,1,1);
%
if max(x(:,4))>0 && max(x(:,4))<45 % AP van
    plot(t, x(:,4),'r')
    drawnow;
    hold on;
    

elseif max(x(:,4))>45 % t�l magas membr�npotenci�l (nem igazi AP)
    plot(t, x(:,4),'b')
    drawnow;
    hold on; %t�bb trajekt�ria egym�sra rajzol�sa
    
else
    plot(t, x(:,4),'k') %nem alakul ki AP
    drawnow;
    hold on;
    
end

    title(['Threshold Behavior, voltage: ', num2str(-65+impulzus)]);
    xlabel('Time (ms)');
    ylabel('Transmembrane Voltage (mV)');
    xlim([-1 50]);

subplot(2,1,2);
%ioncsatorna aktivit�sok alakul�sa
plot(t, x(:,3),'r')
hold on
plot(t, x(:,2),'k')
plot(t, x(:,1),'b')
legend('h: Na gate inactivity','m: Na gate activity','n: K gate activity');
title('Ion gate activities');
xlabel('Time (ms)');
ylabel('Activity');
ylim([0 1]);
xlim([-1 50]);
hold off;

waitforbuttonpress;
impulzus=impulzus+0.2; %az ingerl� potenci�l v�ltoztat�sa iter�ci�nk�nt

end