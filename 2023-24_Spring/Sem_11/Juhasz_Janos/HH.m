% 11. gyakorlat:
% Hodgkin-Huxley modell akciós potenciál kialakulása (1 pont vizsgálata az axon membránon) (ODE rendszer)
% kezdeti (nyugalmi) értékek: ioncsatorna kapuváltozók (mennyire nyitott ~ ellenállás) :
% [K csatorna nyitási, Na csatorna nyitási, Na csatorna inaktiválódási, kezdeti membránpotenciál]

%Feladatok:
% Implementáld a HH egyenletrendszert a hhode.m fájlba!
% Állítsd be a kezdeti ingert (impulzus) (és szükség esetén a többi benenetet), úgy hogy
% A) akciós potenciál (AP) alakuljon ki,
% B) ne legyen AP (marad a nyugalmi potenciál az input után), hiperpolarizált a membrán
% C) ne legyen AP, mert túl magas már a membránpotenciál,
% D) "csípd el" az AP kialakulását (azt az ingerlõ potenciált, ami már épp elég az AP-hez),
% Mi jellemzõ az AP alakjára, megjelenésére, az ioncsatornák aktivitására?
% Mik a nyugalmi értékeik?

figure;
impulzus=; %start impulzus
for i=1:20
ySS=[0.24 0.03 0.76 -65+impulzus]; %nyugalmi potenciál, innen indul a rendszer
[t,x]=ode45('hhode',[0 50],ySS); %%HH egyenlet (külön fájl: hhode.m)
subplot(2,1,1);
%
if max(x(:,4))>0 && max(x(:,4))<45 % AP van
    plot(t, x(:,4),'r')
    drawnow;
    hold on;
    

elseif max(x(:,4))>45 % túl magas membránpotenciál (nem igazi AP)
    plot(t, x(:,4),'b')
    drawnow;
    hold on; %több trajektória egymásra rajzolása
    
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
%ioncsatorna aktivitások alakulása
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
impulzus=impulzus+0.2; %az ingerlõ potenciál változtatása iterációnként

end