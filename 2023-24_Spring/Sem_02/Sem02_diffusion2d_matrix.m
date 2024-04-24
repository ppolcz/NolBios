% 2. gyakorlat: diff�zi� 2 dimenzi�ban (teh�t egy s�kon) m�trixos alak
tic
dt = 0.1;
dx = 1;
dy = 1;

x = 1:dx:100; %1000;
y = 1:dy:100; %1000;
nx = length(x);
ny = length(y);
nt = 1500;

Dx=1;
Dy=1;

u = zeros(nx,ny);
u(70:80, 70:80) = 1; %kezdeti �llapot (n�gyzet)
u(10:20, 70:80) = 1;

figure

for t = 1:nt
    
    %periodikus perem:
    i_before=[u(end, :); u(1:end-1, :)];
    i_after=[u(2:end, :); u(1, :)];
    j_before=[u(:, end), u(:, 1:end-1)];
    j_after=[u(:, 2:end), u(:, 1)];
    
    
    
    u = u + Dx*dt/dx^2*(i_after-2*u+i_before)...
                       + Dy*dt/dy^2*(j_after-2*u+j_before); %dU/dt=D*dU^2/dx^2 %diff�zi�, sz�tfolyik
    
    %imagesc(u);
    mesh(u);
    zlim([0 1.5]); ylim([0 ny]); xlim([0 nx]);
     
    caxis([0 1]); %lecseng�s norm�lt sz�nekkel jobban l�tsz�dik
    %colorbar;
    title(['l�p�s: ', num2str(t)]);
    %waitforbuttonpress;
    pause(0.01);
end
toc