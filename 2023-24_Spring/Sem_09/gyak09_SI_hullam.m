%9. gyakorlat
% SI rendszer egyenlet utaz�hull�ma:
dt = 0.1; %id�beli l�p�sk�z 
dx = 1; %t�rbeli l�p�sk�z

x = 1:dx:1000; %viszg�lt pontok az egyenesen
nx = length(x);
tmax = 150;
nt = round(tmax/dt); %id�tartom�ny v�ge
indv = zeros(nt,1); %sebess�ghez adatok ide;

c = 1; %hul�mterjed�s sebess�gei �lland�
D = 4; %diff�zi� sebess�ge
tau = 2; %fert�z�si r�ta

I = zeros(nx,1); %a betegek kezdeti �rt�kei
%kezdeti �rt�kek:
I(1) = 0.1; %kezdeti fert�z�ttek a t�r elej�n
I(250:300) = 0.001; %1;%%0.5; %0.1; %0.01;% % %kis ter�letr�l indul %fert�z�ttek kezdeti �rt�kt�l nem f�gg a hull�mfront 
%eta(1:100) = linspace(0.3,0,100); %a kezdeti hul�mfrontt�l se f�gg a hull�m alak
S = 1-I; %fert�zhet�k kezdeti �rt�kei
%fi_kezd = fi(1);

figure;
%kezdeti �rt�k rajzol:
plot(1:nx, I, 'r', 1:nx, S, 'g', 'LineWidth', 2);
grid(gca, 'minor');
    ylim([-0.2 1.2]);
    xlim([1 1000]);
    title('piros: fert�z�tt (\eta), z�ld: fert�zhet� (\phi); t=0');
    xlabel('x: 1D t�r');
    waitforbuttonpress;
    
for k = 1:nt %minden id�pillanatban
    S_kp1 = S;
    I_kp1 = I;
    for i = 1:nx %minden pontban kisz�molom a hull�m aktu�lis �rt�k�t (a sz�leket k�l�n kell kezelni)
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
        % dfi/dt=c*fi*(-tau*eta)+D*dfi^2/dx^2 %fert�zhet�ek (S) v�ltoz�sa, reakci�(-beteged�s)+diff�zi�  
        % deta/dt=c*eta*(tau*fi-1)+D*deta^2/dx^2 %fert�z�ttek (I) v�ltoz�sa, reakci�(+beteged�s-gy�gyul�s)+diff�zi� 
        S_kp1(i) = S(i) +c*dt*S(i)*(-tau*I(i))+D*dt/dx^2*(S(utan)-2*S(i)+S(elott));
        I_kp1(i) = I(i) +c*dt*I(i)*(tau*S(i)-1)+D*dt/dx^2*(I(utan)-2*I(i)+I(elott));
    end
    S = S_kp1;
    I = I_kp1;
    
    %�br�zol:
    plot(1:nx, I, 'r', 1:nx, S, 'g', 'LineWidth', 2);
    grid(gca, 'minor');
    ylim([-0.2 1.2]);
    xlim([1 1000]);
    title(['piros: fert�z�tt (\eta), z�ld: fert�zhet� (\phi); t=', num2str(k*dt),]);
    xlabel('x: 1D t�r');
    %waitforbuttonpress;
     pause(0.01);
     
     %{
    % hull�m sebess�g�t meghat�roz: 
 kul=abs(fi-0.5); 
 [mi, hol]=min(kul); %megkeresem a hull�m k�zep�t (ami legk�zelebb van 0.5-h�z)
 indv(t,1)=hol(1); %t�rolom, hogy hol, milyen x-n�l van ez adott id�pillanatban
    %}
end
 
%{
%sebess�g:
elvart_min_sebesseg=2*sqrt(tau*fi_kezd-1) %u(1)=1; t�bbi 0 kezdeti esetben
 tkezd=1/4*nt; %egy intervallumon bel�l, ahol stabil a hull�m sebess�ge
 tveg=4/5*nt;
sebesseg=abs(indv(tveg,1)-indv(tkezd,1))/abs((tveg-tkezd)*dt) %mennyit halad a hull�m k�zepe adott id� alatt 
%}




