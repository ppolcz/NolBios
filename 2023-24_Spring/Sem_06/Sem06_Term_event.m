
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

XLim = [-0.5,5];
YLim = [-0.5,6];
epsilon = 0.01;

% opts = odeset('Events',@termevent);
opts = odeset('Events',@(t,x) termevent_param(t,x,[XLim YLim],Eq,epsilon));

figure,
hold on

r = 0.05;
for eq = Eq
    
    A = J_fh(eq);
    [S,lambda] = eig(A,'vector');

    directions = [S,-S];
    for dir = directions
        x0 = eq + dir * r;
        [t,x] = ode45(f_fh,[0,100],x0,opts);
        plot(x(:,1),x(:,2),'r');
        [t,x] = ode45(f_fh,[0,-100],x0,opts);
        plot(x(:,1),x(:,2),'r');
    
        % plot(x0(1),x0(2),'.b','MarkerSize',5)
    end

end

function [q,isTerm,dir] = termevent_param(t,x,lims,Eq,epsilon)
    % lims = [xmin,xmax,ymin,ymax]

    q = [
        x(1)-lims(1)  % q1
        x(1)-lims(2)  % q2
        x(2)-lims(3)  % q3
        x(2)-lims(4)  % q4
        vecnorm(Eq - x)' - epsilon
        ];
    
    isTerm = [1;1;1;1;ones(size(Eq,2),1)];

    dir = [0;0;0;0;zeros(size(Eq,2),1)];
end

function [q,isTerm,dir] = termevent(t,x)

    q = [
        x(1)+0.5  % q1
        x(1)-5    % q2
        x(2)+0.5  % q3
        x(2)-6    % q4
        ];
    
    isTerm = [1;1;1;1];

    dir = [0;0;0;0];
end