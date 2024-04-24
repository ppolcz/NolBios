
eta_c = 730 \ (26771 + sqrt(80211121));
nu_c = 117 - 2*eta_c;

eta = eta_c + 1/1000;
nu = nu_c + 1/100000;

A = -[
    23  9 eta
    37 16  nu
    9   7  16
    ];
 
b = -sum(A,2);

x = sym('x',[3,1]);
f = diag(x)*(A*x+b);

f_fh = @(t,x) diag(x)*(A*x+b);
J_fh = matlabFunction(jacobian(f,x),'vars',{x});
f1_fh = matlabFunction(f(1),'vars',x);
f2_fh = matlabFunction(f(2),'vars',x);
f3_fh = matlabFunction(f(3),'vars',x);

sol = solve(f,x);
Eq = double([sol.x1,sol.x2,sol.x3]');
Eq(:,any(Eq < 0,1)) = [];

XLim = [-0.5,7];
YLim = [-0.5,7];
ZLim = [-0.5,3];
epsilon = 0.01;

% opts = odeset('Events',@termevent);
opts = odeset('Events',@(t,x) termevent_param(t,x,[XLim YLim ZLim],Eq,epsilon));

figure,
hold on

plot3(Eq(1,:),Eq(2,:),Eq(3,:),'.k','MarkerSize',25)

r = 0.05;
for eq = Eq
    
    A = J_fh(eq);
    [S,lambda] = eig(A,'vector');

    if any(abs(imag(lambda)) > 0)
        S(:,abs(imag(lambda)) > 0) = [];
    end

    directions = [S,-S];
    for dir = directions
        x0 = eq + dir * r;
        [t,x] = ode45(f_fh,[0,100],x0,opts);
        plot3(x(:,1),x(:,2),x(:,3),'r');
        [t,x] = ode45(f_fh,[0,-100],x0,opts);
        plot3(x(:,1),x(:,2),x(:,3),'r');
    
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
        x(3)-lims(5)  % q3
        x(3)-lims(6)  % q4
        vecnorm(Eq - x)' - epsilon
        ];
    
    isTerm = [1;1;1;1;1;1;ones(size(Eq,2),1)];

    dir = [0;0;0;0;0;0;zeros(size(Eq,2),1)];
end
