% 11. gyakorlat:
% Hodgkin-Huxley: akciós potenciál terjedés az axonon (parciális egyenlet is):

% Feladatok:
% Implementáld a parciális egyenletrendszert a hhode_parc.m fájlba!
% Állítsd be a nyugalmi állapot értékeit (elõzõ feladat alapján) és adj ingert ad axon elejére, hogy elinduljon az AP!
% Figyeld meg a terjedõ akciós potenciált!
% A hhode_parc.m fájlban adj különbözõ kezdõpontokban különbözõ méretû 2. impulzust is az axon elejére, és figyeld meg, mikor alakul ki AP!

close all;
% kezdeti (nyugalmi) értékek: ioncsatorna kapuváltozók (mennyire nyitott ~ ellenállás) :
% [K csatorna nyitási, Na csatorna nyitási, Na csatorna inaktiválódási, kezdeti membránpotenciál]
resting = [1 2 4 5]; % nyugalmi értékek
firstpot = 1; % kezdeti potenciálváltozás axon elején

N = 60; % az axon hossza (ennyi egységbõl áll)
% az axon minden pontjához a vektor 4 mezõjét rendelem (n, m, h, V), 
% tehát egy pontot egy 4 hosszú vektor jelöl, 
% ezért léptetek az ySS mentén 4-esével
ySS = zeros(1,N*4);
% kezdetben az axon minden pontja nyugalomban van 
ySS(1:4:end) = resting(1);
ySS(2:4:end) = resting(2);
ySS(3:4:end) = resting(3);
ySS(4:4:end) = resting(4);
ySS(4) = firstpot; %az 1. pozíció membránpotenciálja megváltozik -> ingerület 

[t,x]=ode15s('hhode_parc',[0,60],ySS); %%HH egyenlet (külön fájl: hhode_parc.m) 
% ose15s-sel stabilabban mûködik, (stiff probléma a hiperpolarizáció alatti gerjesztés...?)

%
% megjelenítés: 
figure;
for i = 1:length(t)
   %plot(1:N, x(i,4:4:4*N));
   plot(linspace(1,N,400),  interp1(x(i,4:4:4*N),linspace(1,N,400),'pchip')  ); %a hhode_parc-ban kiszámolt értékeket (x minden 4. eleme a membránpotenciál) finomítom a linspace-szel megadott pontossághoz
  % 400 helyett pl 20 ponton interpolálok-> darabosabb fv, a pulzálás a numerikus pontatlanság miatt van
  xlim([0, N]);
  ylim([-90, 60]);
  title(['Time: ', num2str(t(i))]);
  xlabel('Axon');
  ylabel('Activity');
  pause(0.005);
% waitforbuttonpress;
end
%}