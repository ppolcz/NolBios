

syms t c x1 x2 x3 real
x = [x1;x2;x3];

a = 0.1;
b = 0.1;
% c = 5.7; % from 4 to 18

f = [
    -x2 - x3
    x1 + a*x2
    b + x3*(x1-c)
    ];
f_ode = matlabFunction(f,'vars',{t,x,c});

%%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

x0 = [0.1;3;0.11];
T_min = 100;
T = 1000;

fig = figure(123);
fig.Position(3:4) = [1421 758];
Tl = tiledlayout(3,2,'Padding','tight');
ax = nexttile([3,1]); hold on, grid on, box on;
Pl = plot3(0,0,0);

view([5.6 13]);
xlabel x1
ylabel x2
zlabel x3
axis equal

ax.XLim = [-30 30];
ax.YLim = [-30 30];
ax.ZLim = [0 60];

clim = [4,18];

nexttile, hold on, grid on, box on;
PlX = plot(0,0,'.','MarkerSize',2);
xlim(clim);
ylim([0,30])
xlabel('Peaks in X')

nexttile, hold on, grid on, box on;
PlY = plot(0,0,'.','MarkerSize',2);
xlim(clim);
ylim([0,25])
xlabel('Peaks in Y')

nexttile, hold on, grid on, box on;
PlZ = plot(0,0,'.','MarkerSize',2);
xlim(clim);
ylim([0,60])
xlabel('Peaks in Z')


%%

for c_val = linspace(clim(1),clim(2),500)
    c_val
    [tt,xx] = ode89(@(t,x) f_ode(t,x,c_val),[0,T],x0);
    xx = xx(tt > T_min,:);

    Pl.XData = xx(:,1);
    Pl.YData = xx(:,2);
    Pl.ZData = xx(:,3);

    [PEAKS,LOCS] = findpeaks(xx(:,1));
    PlX.XData = [PlX.XData , PEAKS'*0+c_val];
    PlX.YData = [PlX.YData , PEAKS'];

    [PEAKS,LOCS] = findpeaks(xx(:,2));
    PlY.XData = [PlY.XData , PEAKS'*0+c_val];
    PlY.YData = [PlY.YData , PEAKS'];

    [PEAKS,LOCS] = findpeaks(xx(:,3));
    PlZ.XData = [PlZ.XData , PEAKS'*0+c_val];
    PlZ.YData = [PlZ.YData , PEAKS'];

    drawnow
end


