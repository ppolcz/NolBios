%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. November 09. (2023a)
%  Revised on 2025. October 07. (2025a)

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

fig = figure(2);
delete(fig.Children)
fig.Color = [1 1 1];
Tl = tiledlayout(2,2,"TileSpacing","compact","Padding","compact","TileIndexing","rowmajor");

Ax_Bif = nexttile(2);
hold on; grid on; box on;
Tit_Bif = title("$\mu = " + 0 + "$",'Interpreter','latex','FontSize',15);
xlabel('$\mu$ (the parameter)','Interpreter','latex','FontSize',12)
ylabel('$x_{\infty}$ (the "steady states")','Interpreter','latex','FontSize',12)

X = [];
Y = [];

for mu = linspace(1,4,301)
    f = @(x) mu*(x - x.^2);
    
    N0 = 1000;
    N = 5000;
    x = [ 0.5 zeros(1,N-1) ];
    
    for k = 1:N-1
        x(k+1) = f(x(k));
    end
    
    w = limit_points(x(N0:end));

    X = [X ; w(:)*0+mu];
    Y = [Y ; w(:)];
end

Pl_Bif = plot(X,Y,'.');
Pl_mu = plot([0,0],[0,1],'r','LineWidth',2);

Ax_Log = nexttile(1);
hold on; grid on; box on;
Tit_Web = title('$x_{n+1} = \mu \, (x_n - x_n^2)$, where $mu = 0$','Interpreter','latex','FontSize',15);
xlabel('$x_n$','Interpreter','latex','FontSize',12)
ylabel('$x_{n+1}$','Interpreter','latex','FontSize',12)

t = 0:0.01:1;
plot(t,t,'Color',Color_5,'LineWidth',1.5,'DisplayName','$y = x$');
Pl_f = plot(t,f(t),'Color',Color_1,'LineWidth',1.5,'DisplayName','$y = \mu\,(x - x^2)$');

[X,Y] = plotXY(x,f);
Pl_WEB = plot(X,Y,'Color',Color_1,'HandleVisibility','off');

[X,Y] = plotXY(x(N0:end),f);
Pl_W = plot(X(2:end),Y(2:end),'Color',Color_2,'LineWidth',2,'DisplayName','$x_\infty$ (steady state)');

Leg = legend('Interpreter','latex','FontSize',12,'Location','northwest');

Felosztas = 0:0.01:1;
Ax_Hist = nexttile(3);
hold on; grid on; box on;
Pl_Hist = histogram(x(N0:end),Felosztas);
title("Distribution of ""steady state"" values",'Interpreter','latex','FontSize',15);
xlabel('$x_n$','Interpreter','latex','FontSize',12)
ylabel('occurrence','Interpreter','latex','FontSize',12)
xlim([0,1])

Ax_PLOT = nexttile;
hold on; grid on; box on;
Pl_idofv = stairs(1:N,x);
xlim([0,50])
title("Trajectory $x_n$, where $n = 0,1,2,\dots$",'Interpreter','latex','FontSize',15);
xlabel('time ($n$)','Interpreter','latex','FontSize',12)
ylabel('$x_n$','Interpreter','latex','FontSize',12)
ylim([0,1]);

%%

tic
for mu = [ linspace(0,1,51) linspace(1,3.5,301) linspace(3.5,4,301) ]
    Tit_Bif.String = "$\mu = " + mu + "$";
    Tit_Web.String = "$x_{n+1} = \mu \, (x_n - x_n^2)$, where $\mu =" + mu +"$";
    
    f = @(x) mu*(x - x.^2);
    
    N0 = 1000;
    N = 5000;
    x = [ 0.01 zeros(1,N-1) ];
    
    for k = 1:N-1
        x(k+1) = f(x(k));
    end
    w = limit_points(x(N0:end));

    Pl_f.YData = f(t);
    [Pl_WEB.XData,Pl_WEB.YData] = plotXY(x,f);
    [Pl_W.XData,Pl_W.YData] = plotXY(x(N0:end),f);
    
    Pl_mu.XData = mu*[1,1];
    % Pl_Hist = histogram(Ax_Hist,x(N0:end),Felosztas);
    Pl_Hist.Data = x(N0:end);

    Pl_idofv.YData = x;

    pause(0.1-toc)
    tic

    if ismember(mu,[2,3,3.5,3.7,3.83,3.9])
        keyboard
    end
end

function [X,Y] = plotXY(x,f)
    N = numel(x);
    X = reshape([1;1]*x, [2*N 1]);
    Y = reshape([ 0 f(x(1:end-1)) ; f(x) ],[2*N 1]);
end

function w = limit_points(x)
    w = unique(round(x,10));
end