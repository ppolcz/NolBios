

syms t rho x1 x2 x3 real
x = [x1;x2;x3];

% rho = 28;
sigma = 10;
beta = 8/3;

f = [
    sigma*(x2 - x1)
    x1*(rho - x3) - x2
    x1*x2 - beta*x3
    ];

f_ode = matlabFunction(f,'vars',{t,x,rho});

%%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

x0 = [0.1;3;0.11];
T = 100;

fig = figure(123);
delete(fig.Children)
ax = axes(fig); hold on, grid on, box on;
Pl = plot3(0,0,0);

view([5.6 13]);
xlabel x1
ylabel x2
zlabel x3
axis equal

%%

for rho_val = linspace(20,28,1000)
    rho_val
    [tt,xx] = ode89(@(t,x) f_ode(t,x,rho_val),[0,T],x0);

    Pl.XData = xx(:,1);
    Pl.YData = xx(:,2);
    Pl.ZData = xx(:,3);

    drawnow
end


