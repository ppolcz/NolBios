%9. gyakorlat:
% Fisher egyenlet 2 elsõrendû egyenletre bontva, fázisportré ábrázolva
% utazó hullám egyenlet: -s*dx=x*(1-x)+ddx ->
% dx=y
% dy=-s*y-x*(1-x)=ddx

figure;
x=-0.1:0.01:1.1;% -1.1:0.01:2.1;%
y=-0.3:0.01:0.2;% -1.3:0.01:1.2;%
[X,Y]=meshgrid(x,y);

s=1;%3;%1;%2; %hullám sebességi állandó % 2 alatt fókusz van 0-ban, biológiailag nem releváns, mert negatív valséget jelentene
% vektormezõs ábrázolás
dx=Y; 
dy=-s*Y-X.*(1-X);
%l=sqrt(dx.^2+dy.^2); %normálás
%dx=dx./l ; 
%dy=dy./l; 
quiver(X,Y,dx, dy, 2, 'k'); 
hold on;

%egy pálya kirajzolása:
%{
Xinit=0.9;
Yinit=-0.01;
tint=[0 100];
hullam_ode=@(t,xy)[xy(2); -s*xy(2)-xy(1)*(1-xy(1))];
[t, xy] = ode45(hullam_ode,tint,[Xinit Yinit]); %az diffegyenletet kiszámolom
plot(xy(:,1),xy(:,2),'b');
%}

%több pálya kirajzolása:
%
N=50;
%az ábra 4 határoló egyeneseirõl indítok trajektóriákat:
sxk=x(1);
syk=y(1);
sxv=x(1);
syv=y(end);

sx2k=x(1);
sy2k=y(1);
sx2v=x(end);
sy2v=y(1);

sx3k=x(end);
sy3k=y(1);
sx3v=x(end);
sy3v=y(end);

sx4k=x(1);
sy4k=y(end);
sx4v=x(end);
sy4v=y(end);

sx=[linspace(sxk,sxv,N), linspace(sx2k,sx2v,N), linspace(sx3k,sx3v,N), linspace(sx4k,sx4v,N)];
sy=[linspace(syk,syv,N), linspace(sy2k,sy2v,N), linspace(sy3k,sy3v,N), linspace(sy4k,sy4v,N)];

%adott pontokból induló trajeltóriák kiszámítása:
tint=[0,100];
for i=1:length(sx)
% aktuális kezdõpontok
Xinit=sx(i);
Yinit=sy(i);
hullam_ode=@(t,xy)[xy(2); -s*xy(2)-xy(1)*(1-xy(1))]; %differenciál egyenlet
[t, xy] = ode45(hullam_ode,tint,[Xinit Yinit]); %a diffegyenletet kiszámolom
plot(xy(:,1),xy(:,2),'r');
end
%}

xlabel('x');
ylabel('y');
xlim([x(1)-0.1;x(end)+0.1]);
ylim([y(1)-0.1;y(end)+0.1]);
