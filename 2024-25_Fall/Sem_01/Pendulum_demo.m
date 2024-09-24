%%
% Author: Peter Polcz (ppolcz@gmail.com) 
% Created on 2023. September 28. (2023a)

%% Damped oscillation, no external input

m = 1;
l = 1;
g = 10;
b = 0.2;

% function ret = pendulum_ODE(t,x)
%     ret = [
%         x(2)
%         -g/l*sin(x(1)) - b/m*x(2)
%         ];
% end
% 
% This is actually a function, declared in a simplified form.
% 
Eq = @(t,x) [
    x(2)
    -g/l*sin(x(1)) - b/m*x(2)
    ];

% This function can be simply called as follows:
t = 13;
x = [2;3];
Eq(t,x)

% Energy function
V_fh = @(x1,x2) m*g*l*(1-cos(x1)) + 1/2*m*l^2*x2.^2;

% This function can be simply called as follows:
x1 = 12;
x2 = 1;
V_fh(x1,x2)

x0 = [
    -3*pi  % = x1(0)
     6     % = x2(0)
    ];
T = 20;    % The solution is computed for t = [0,T]
[t_sol,x_sol] = ode45(Eq,[0,T],x0);

hp_animate_pendulum(t_sol,x_sol,V_fh,'V_Levels',m*g*l * (0:0.5:5),...
    'Record',false,'ExportName','Pendulum_damped_oscillation')

%% Damped oscillation, periodic external input

m = 1;
l = 1;
g = 10;
b = 0.03;

A = 4;
w = 1.8;
% w = 50;
u = @(t) A*sin(w*t);

Eq = @(t,x) [
    x(2)
    -g/l*sin(x(1)) - b/m*x(2) + u(t)
    ];

% Eq = @(t,x) [
%     x(2)
%     -g/l*x(1) - b/m*x(2) + u(t)
%     ];

% Energy function
V_fh = @(x1,x2) m*g*l*(1-cos(x1)) + 1/2*m*l^2*x2.^2;

x0 = [-3*pi;0];
T = 100;
opts = odeset('MaxStep',0.02);
opts = odeset();
[t_sol,x_sol] = ode45(Eq,[0,T],x0,opts);

% t_eqd = linspace(0,T,numel(t_sol));
% y = interp1(t_sol,x_sol(:,1),t_eqd);
% Y = abs(fft(y,2^13));
% U = abs(fft(u(t_eqd),2^13));
% plot(U), hold on
% plot(Y)

hp_animate_pendulum(t_sol,x_sol,V_fh,u,'V_Levels',m*g*l * (0:0.5:5),...
    'TeXInput',sprintf('%g\\sin(%g t)',A,w),...
    'Record',false,'ExportName','Pendulum_chaotic_behaviour')

%%

% Bolzano shooting

% Parameters:
X0 = 0;         % Initial conditions
Y0_u = 4.9;     % Upper
Y0_l = 1;       % Lower
t0 = 0;         % start time
T = 100;        % end time point
b = 0.2;        % damping parameter
eps = 1.0e-6;   % required accuracy
goal_x = -pi;   % target value (get close to this)
goal_y = 0;

%make grid, calculate the energy levels
x = linspace(-pi,3*pi,1000);
y = linspace(-3,Y0_u,500);
[X,Y] = meshgrid(x,y);
Z = (Y.^2)/2 + (1 - cos(X)); 

% Equation of damped pendulum
Eq = @(t,x) [
     x(2)
    -b*x(2)-sin(x(1))
    ];

V_fh = @(x1,x2) (1-cos(x1)) + 1/2*x2.^2;

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];
Xis2Pi_Color = [0.8 0.1 0.1];
Xis0_Color = [0.3 0.1 0.5];

fig = figure(1236);
fig.Position(3:4) = [1159,617];
Tl = tiledlayout(1,2,"TileSpacing","compact","Padding","tight");
Ax = nexttile; 
hold on
xlim(x([1,end]))
ylim(y([1,end]))

plot([X0 X0],[Y0_l Y0_u],'ko-')
contour(X,Y,Z)

for It = 1:1000 
    % midpoint:
    Y0 = (Y0_u + Y0_l)/2; 
    
    % solve the equation
    [t,x_sol] = ode45(Eq,[t0 T],[X0;Y0]); 

    if x_sol(end,1) <= pi

        Y0_l = Y0;
        Color = [0 0 1];

    else

        Y0_u = Y0; 
        Color = [1 0 0];

    end

    % Draw the first few result
    if It < 10
        plot(x_sol(1,1),x_sol(1,2),'.','Color',Color,'MarkerSize',15)
        drawnow
        
        pause(0.1)
    
        plot(x_sol(:,1),x_sol(:,2),'Color',Color,'LineWidth',1.5);
        drawnow
    
        pause(0.1)
    end
end

% Draw the last result
plot(x_sol(1,1),x_sol(1,2),'.','Color',Color,'MarkerSize',15)
plot(x_sol(:,1),x_sol(:,2),'Color',Color,'LineWidth',1.5);
drawnow

Ax = [Ax nexttile]; 
hold on

plot([X0 X0],[Y0_l Y0_u],'ko-')
contour(X,Y,Z)

[~,Idx] = min(abs(x_sol(:,1) - pi));
x_sol = x_sol(1:Idx,:);

t_uniform = t0:0.1:6; %t(Idx);
x_sol = interp1(t(1:Idx),x_sol,t_uniform);
plot(x_sol(:,1),x_sol(:,2),'.')

xlim(x([1,end]))
ylim(y([1,end]))

for ax = Ax
    grid(ax,'on'); 
    box(ax,'on');
    
    Font = {'Interpreter','latex','FontSize',16};
    xlabel('$x$ [rad]',Font{:});
    ylabel('$y$ [rad/s]',Font{:});
    
    ax.TickLabelInterpreter = "latex";
    ax.FontSize = Font{4};
    
    x_grid_pi_mtp = -10:10;
    XTickLabels = cellfun(@(s) {"$" + strrep(latex(s),'\,','') + "$"},num2cell(x_grid_pi_mtp*sym(pi)));
    ax.XTick = x_grid_pi_mtp * pi;
    ax.XTickLabel = XTickLabels;
end

ax = nexttile('south');
plot(t_uniform,x_sol);
hold on, grid on
Leg = legend('Position ($x_1$)','Velocity ($x_2$)',...
    'Interpreter','latex','FontSize',14);
title(sprintf('The velocity needed to flip the pendulum up from is $%g$',x_sol(1,2)),...
    'Interpreter','latex','FontSize',14)

exportgraphics(fig,'Pendulum_Compute_Swingup_speed_Bolzano.png','ContentType','image');

hp_animate_pendulum(t_uniform,x_sol,V_fh,'X2Lims',[-3,3],...
    'Record',false,'ExportName','Pendulum_swingup')
