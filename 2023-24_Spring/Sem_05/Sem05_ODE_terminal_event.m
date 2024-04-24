%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. April 10. (2023a)
%

syms t x y mu eta real

mu = 1.5;
eta = 2;

f = [
    2*x*(1 - x) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

sol = solve(f);
Eq = double([
    sol.x sol.y
    ]');

f_fh = matlabFunction(f,'vars',{t,x});           % later used for ODE solver
J_fh = matlabFunction(jacobian(f,x),'vars',{x}); % compute trace and determinant of Jacobian
f1_fh = matlabFunction(f(1),'vars',x);           % to plot vector field
f2_fh = matlabFunction(f(2),'vars',x);           % to plot vector field

str = "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";

XLim = [-0.5,4];
YLim = [-0.5,4];

termevent = @(t,x) hp_ode_terminal_event_pp(t,x,XLim,YLim,Eq,0.01);
opts = odeset('Events',termevent);

fig = figure(12);
delete(fig.Children)
ax = axes(fig);
hold on, grid on, box on;

% Vektormezo
[x,y] = meshgrid(linspace(XLim(1),XLim(2),101),linspace(YLim(1),YLim(2),101));

f1 = f1_fh(x,y);
f2 = f2_fh(x,y);

r = sqrt(f1.^2 + f2.^2);
f1 = f1 ./ r;
f2 = f2 ./ r;

% quiver(x,y,f1,f2)
streamslice(x,y,f1,f2)


r = 0.02;
for eq = Eq

    [S,lambda] = eig(J_fh(eq),'vector');

    for x0 = eq + [S , -S] * r

        [t,x] = ode45(f_fh,[0,100],x0,opts);
        plot(x(:,1),x(:,2),'r','LineWidth',2)
        [t,x] = ode45(f_fh,[0,-100],x0,opts);
        plot(x(:,1),x(:,2),'r','LineWidth',2)
    end
end

xlim(XLim);
ylim(YLim)