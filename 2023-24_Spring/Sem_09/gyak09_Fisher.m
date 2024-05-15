%9. gyakorlat:
% Fisher egyenlet:
dt = 0.1; %id�beli l�p�sk�z %D*2*dt/dx^2<1 %0.46 k�r�l "sz�p sz�r�s�d�s", nagy ,kis freki ugyan�gy cseng le
dx = 1; %t�rbeli l�p�sk�z

x = 1:dx:500; %viszg�lt pontok az egyenesen
nx = length(x);
tmax=150;%1500;%75;%
nt = round(tmax/dt); %id�tartom�ny v�ge
indv=zeros(nt,1); %sebess�ghez adatok ide;

c =1; %hul�mterjed�s sebess�gi �lland�ja
D=1; %diff�zi� sebess�ge

u=zeros(nx,1);%exp(-0.1*x);%  %a hull�m kezdeti �rt�kei (ha sehol se 0 kezdetben a rendszer, az befoly�solja a hull�mfrontot) 
%kezdeti �rt�kek:
u(1) =1; %5;%0.1;%1;% az 1. pontban 1 (a t�bbiben 0) a hull�m kezdetben
u(end) = 1;
% u(1:50) =1;%0.1;%5;%1; % (pozit�v) kezdett�l f�ggetlen�l 1-re �ll be a hull�m
% u(250:end)=1;% visszafele is terjed
% u(250:260)=5;%5;%0.1;%1; %cs�csb�l indul 
% u([250:260, 200:230])=1; %0.1;%5;%1;% 2 cs�cs indul 
% u(250:300)=rand(1,51)*2; %random kezdet (nem lehet negat�v)
% u(250:300)=linspace(0,2,51); %a kezdeti hul�mfrontt�l se f�gg hull�m alakja
% u(1:100)=linspace(2,0,100); %a kezdeti hul�mfrontt�l se f�gg
% u(1:100)=linspace(0.5,0,100); %a kezdeti hul�mfrontt�l se f�gg
% n�gyzetets Fisher:
% u(250:300)=rand(1,51)*2-1; %random kezdet (lehet negat�v)
u(250:300)=linspace(-0.5,2,51); %a kezdeti hul�mfrontt�l se f�gg (negat�v r�sz is)
%u(200:249)=linspace(-1,1,50); u(250:300)=linspace(1,-0.1,51); u(301:350)=linspace(-0.1,1,50);%majdnem -1-re �ll be
%u(200:249)=linspace(-1,1,50); u(250:300)=linspace(1,-0.06,51); u(301:350)=linspace(-0.06,1,50);%sok�ig �gy t�nik -1-re �ll be, azt�n ugrik fel 1-re (hosszabb fut�s): metastabilit�s 

figure;
%kezdeti �rt�k kirajzol�sa:
plot(u); %hull�m �rt�keinek �br�zol�sa
grid(gca, 'minor');
    ylim([-1.2 5]);
    xlim([1 450]);
    title('t=0');
    xlabel('x: 1D t�r');
    waitforbuttonpress;
    
for t = 1:nt %minden id�pillanatban
    u2 = u; %u2=�rt�kek az el�z� id�plillanatban; u=�j �rt�kek
    for i = 1:nx %minden pontban kisz�molom a hull�m aktu�lis �rt�k�t (a sz�leket k�l�n kell kezelni)
        %kezelj�k a peremeket zero flux peremk�nt
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
        % 1. dU/dt = c*U*(1-U) + D*dU^2/dx^2 % reakc tag+diff tag=fisher, t�rben sz�tterjed (diff),
        % u(i) = u2(i) + c*dt*u2(i)*(1-u2(i)) + D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % 2. dU/dt = c*U*(1-U^2) + D*dU^2/dx^2 % n�gyzetes fisher (-1-es stabil pont is van)
        u(i) = u2(i) +c*dt*(u2(i))*(1-u2(i)^2)+D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % hasonl� rendszerek:
        % 3. dU/dt=D*dU^2/dx^2 %diff tag csak, sz�tfolyik ,de nem mozdul
        % u(i) = u2(i) + D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        % 4. dU/dt=-c*dU/dx %lin convection: t�rben utazik, nem folyik sz�t a hull�m  �s fenn is marad (reakc) 
        % u(i) = u2(i) - c*dt/dx*(u2(i)-u2(elott));
        % 5. dU/dt=-c*dU/dx+D*dU^2/dx^2 %lin convection+diff=burgers, utazik �s sz�tfolyik 
        % u(i) = u2(i) - c*dt/dx*(u2(i)-u2(elott))+D*dt/dx^2*(u2(utan)-2*u2(i)+u2(elott)); 
        
    end
    
    %�br�zol:
    plot(u);
    grid(gca, 'minor');
    ylim([-1.2 5]);
    xlim([1 450]);
    title(['t=' ,num2str(t*dt)]);
    xlabel('x: 1D t�r');
    %waitforbuttonpress;
     pause(0.01);
     
     %
 % fishern�l hull�m sebess�g�t meghat�roz: (1-b�l ind�tott hull�mn�l m�k�dik j�l)
 kul=abs(u-0.5); 
 [mi, hol]=min(kul); %megkeresem a hull�m k�zep�t (ami legk�zelebb van 0.5-h�z)
 indv(t,1)=hol(1); %t�rolom, hogy hol, milyen x-n�l van ez adott id�pillanatban
    %}
end
 
%
%fishern�l:
elvart_sebesseg=2*sqrt(c*D) %u(1)=1; t�bbi 0 kezdeti esetben
tkezd=round(1/4*nt); %egy intervallumon bel�l, ahol stabil a hull�m sebess�ge (,hogy a tranziensek, peremek ne zavarjanak) 
tveg=round(4/5*nt);
sebesseg=(indv(tveg,1)-indv(tkezd,1))/((tveg-tkezd)*dt) %mennyit halad a hull�m k�zepe adott id� alatt 
%nem lesz pont az elv�rt; 2.036-ra 2 helyett kihozhat� c=1 D=1, dx=1, dt=0.01 eset�n
%}




