
% Damping parameter
b = 0.1;      

% Equation of damped pendulum
Eq = @(t,x) [
     x(2)
    -b*x(2)-sin(x(1))
    ];

% Sebesseg kezdeti pontja (a felso hatarra vonatkozoan)
Y0f = linspace(3.1,4.2,500);  
% Pozicio kezdeti pontja (az felse hatarra vonatkozoan)
X0f = Y0f*0 - 6*pi;

% Sebesseg kezdeti pontja (a also hatarra vonatkozoan)
Y0a = linspace(-4,-3.1,500);  
% Pozicio kezdeti pontja (az also hatarra vonatkozoan)
X0a = Y0a*0 + 6*pi;

initial_conditions = [
    [X0a ; Y0a] [X0f ; Y0f] 
    ]; 


fig = figure(321);
delete(fig.Children)
ax = axes(fig);
hold on;
grid on;
box on;

T = 100;
for x0 = initial_conditions

    [t,x_sol] = ode45(Eq,[0 T],x0);

    Color = [0 0 0];
    if -3*pi <= x_sol(end,1) && x_sol(end,1) <= -pi
        Color = [1 0 0];
    elseif -pi < x_sol(end,1) && x_sol(end,1) <= pi
        Color = [0 1 0];
    elseif pi < x_sol(end,1) && x_sol(end,1) <= 3*pi
        Color = [0 0 1];
    end

    plot(x_sol(:,1),x_sol(:,2),'Color',Color)
    drawnow
end



