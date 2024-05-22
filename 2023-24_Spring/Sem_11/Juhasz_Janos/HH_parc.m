% 11. gyakorlat:
% Hodgkin-Huxley: akci�s potenci�l terjed�s az axonon (parci�lis egyenlet is):

% Feladatok:
% Implement�ld a parci�lis egyenletrendszert a hhode_parc.m f�jlba!
% �ll�tsd be a nyugalmi �llapot �rt�keit (el�z� feladat alapj�n) �s adj ingert ad axon elej�re, hogy elinduljon az AP!
% Figyeld meg a terjed� akci�s potenci�lt!
% A hhode_parc.m f�jlban adj k�l�nb�z� kezd�pontokban k�l�nb�z� m�ret� 2. impulzust is az axon elej�re, �s figyeld meg, mikor alakul ki AP!

close all;
% kezdeti (nyugalmi) �rt�kek: ioncsatorna kapuv�ltoz�k (mennyire nyitott ~ ellen�ll�s) :
% [K csatorna nyit�si, Na csatorna nyit�si, Na csatorna inaktiv�l�d�si, kezdeti membr�npotenci�l]
resting = [1 2 4 5]; % nyugalmi �rt�kek
firstpot = 1; % kezdeti potenci�lv�ltoz�s axon elej�n

N = 60; % az axon hossza (ennyi egys�gb�l �ll)
% az axon minden pontj�hoz a vektor 4 mez�j�t rendelem (n, m, h, V), 
% teh�t egy pontot egy 4 hossz� vektor jel�l, 
% ez�rt l�ptetek az ySS ment�n 4-es�vel
ySS = zeros(1,N*4);
% kezdetben az axon minden pontja nyugalomban van 
ySS(1:4:end) = resting(1);
ySS(2:4:end) = resting(2);
ySS(3:4:end) = resting(3);
ySS(4:4:end) = resting(4);
ySS(4) = firstpot; %az 1. poz�ci� membr�npotenci�lja megv�ltozik -> inger�let 

[t,x]=ode15s('hhode_parc',[0,60],ySS); %%HH egyenlet (k�l�n f�jl: hhode_parc.m) 
% ose15s-sel stabilabban m�k�dik, (stiff probl�ma a hiperpolariz�ci� alatti gerjeszt�s...?)

%
% megjelen�t�s: 
figure;
for i = 1:length(t)
   %plot(1:N, x(i,4:4:4*N));
   plot(linspace(1,N,400),  interp1(x(i,4:4:4*N),linspace(1,N,400),'pchip')  ); %a hhode_parc-ban kisz�molt �rt�keket (x minden 4. eleme a membr�npotenci�l) finom�tom a linspace-szel megadott pontoss�ghoz
  % 400 helyett pl 20 ponton interpol�lok-> darabosabb fv, a pulz�l�s a numerikus pontatlans�g miatt van
  xlim([0, N]);
  ylim([-90, 60]);
  title(['Time: ', num2str(t(i))]);
  xlabel('Axon');
  ylabel('Activity');
  pause(0.005);
% waitforbuttonpress;
end
%}