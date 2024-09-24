%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         Practice 1          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonlinear Dynamical Systems %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         2023.09.14.         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

% Parameters:
tinit =0;                   % initial time              % These parameters can be changed
Max_Time = 100;             % end of the simulation     % What happens in longer simulations?
b = 0.1;%0.5;%1.98;%2.02;%0.01;%3;% 2;%           % damping                   % Change the damping value! What happens?
 
x =-3*pi:0.05:5.01*pi;            % not necessary to modify
y =-4:0.05:4;              % not necessary to modify
[X,Y] = meshgrid(x,y);     % makes grid based on x and y vectors (for calculating the energies in each point)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          Exercises          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How do the systems behave? 
% What are the effecs of the parameters? 
%Try other initial condition and plot the into the same plot!
% Are there stable points in the systems? Zoom into the fix points!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          equation 1         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_1 = @(t,x) [
     x(2)
    -x(1)
    ]; 

energy_string = @(X,Y) (Y.^2)/2 + (X.^2)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 2          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_2 = @(t,x) [
     x(2)
    -b*x(2)-x(1)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 3          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_3 = @(t,x) [
     x(2)
    -x(1)+cos(t)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 4          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter (arbitrary), possible values: 0.5, 1.5, etc.
w = 0.5;

eq_4 = @(t,x) [
     x(2)
    -x(1) + cos(w*t)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 5          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_5 = @(t,x) [
     x(2)
    -2*x(1)+cos(t)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 6          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameter (arbitrary), possible values: 0.1, 0.5, 1, 1.5, etc.
w = 1;

eq_6 = @(t,x) [
     x(2)
    -b*x(2)-x(1)+cos(w*t)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 7          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_7 = @(t,x) [
     x(2)
    -sin(x(1))
    ];

energy_pendulum = @(X,Y) (Y.^2)/2 + (1 - cos(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 8          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_8 = @(t,x) [
     x(2)
    -b*x(2)-sin(x(1))
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 9          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eq_9 = @(t,x) [
     x(2)
    -sin(x(1))-b*x(2)+cos(t)
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         equation 10         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = 0.9;
b = 0.05;
eq_10 = @(t,x) [
     x(2)
    -sin(x(1))-b*x(2)+cos(w*t)
    ];

Z = energy_pendulum(X,Y); % Energy (string or pendulum)

Xinit = -5;%-2;                 % initial value of x        %
Yinit = 3;%5;%2;                  % initial value of y        %

% 2. initial point
Xinit2 = -5;%3;%-2.05;%-0.2007;
Yinit2 = 2.99;%2.05;%0.3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             ODE45              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = odeset('RelTol',1e-3, 'AbsTol', 1e-24);

[t,x_sol1]=ode45(eq_9, [tinit Max_Time], [Xinit; Yinit],opts);
[t2,x_sol2]=ode45(eq_9, [tinit Max_Time], [Xinit2; Yinit2],opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         PLOTTING          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        x-y , energy       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fig = figure(1231);
fig.Position(3:4) = [1027,645];
Tl = tiledlayout(3,2,"TileSpacing","compact","Padding","compact","TileIndexing","columnmajor");
Ax1 = nexttile; hold on; grid on; box on;
energy_levels=[(linspace(sqrt(0),sqrt(2),5)).^2, (linspace(sqrt(2),sqrt(8),5)).^2];%linspace(0,8,9);                             % energy levels to calculate and plot
%Choose better energy lines to visualise (linear scale, visualise fix points)!
contour(X,Y,Z, energy_levels);                             % plot chosen energy levels
plot(x_sol1(:,1), x_sol1(:,2), 'b');            
plot(x_sol2(:,1), x_sol2(:,2), 'r');
xlabel('$x$ [rad]');
ylabel('$y$ [rad/s]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             x-t           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Ax2 = nexttile([2,1]); hold on; grid on; box on;
plot(x_sol1(:,1), t, 'b');
plot(x_sol2(:,1), t2, 'r');
xlabel('$x$ [rad]');
ylabel('$t$ [s]');
view([0,-90])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             t-y          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ax3 = nexttile; hold on; grid on; box on;
plot(t, x_sol1(:,2), 'b');
plot(t2, x_sol2(:,2), 'r');
xlabel('$t$ [s]');
ylabel('$y$ [rad/s]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       x-y-t 3D plot     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ax4 = nexttile([2,1]); hold on; grid on; box on;
plot3(x_sol1(:,1), x_sol1(:,2), t, 'b');
plot3(x_sol2(:,1), x_sol2(:,2), t2, 'r');
xlabel('$x$ [rad]');
ylabel('$y$ [rad/s]');
zlabel('$t$ [rad/s]');
Ax4.View = [-15,12];
axis vis3d
rotate3d on
%}

LinkX = linkprop([Ax1,Ax2,Ax4],'XLim');
LinkY = linkprop([Ax1,Ax3,Ax4],'YLim');

for ax = [Ax1,Ax2,Ax3,Ax4]
    ax.TickLabelInterpreter = "latex";
    ax.FontSize = 12;
    ax.XLabel.Interpreter = "latex";
    ax.YLabel.Interpreter = "latex";
    ax.ZLabel.Interpreter = "latex";
end

for ax = [Ax1,Ax2,Ax4]

    x_grid_pi_mtp = -10:10;

    XTickLabels = cellfun(@(s) {"$" + strrep(latex(s),'\,','') + "$"},num2cell(x_grid_pi_mtp*sym(pi)));

    ax.XTick = x_grid_pi_mtp * pi;
    ax.XTickLabel = XTickLabels;

end

