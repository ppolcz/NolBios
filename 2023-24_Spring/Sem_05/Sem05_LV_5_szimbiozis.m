%% Kompetitív Lotka-Volterra modell
% Két egymást segító populáció dinamikájának modellje.
% 

syms t x y mu real
f = [
    1*x*(-1) + x*y
    2*y*(1 - y/2) + x*y
    ];
x = [x;y];
str = "$\dot x = " + latex(f(1)) + "$,~~ $\dot y = " + latex(f(2)) + "$";

A = double(jacobian( simplify(f ./ x) , x ));
b = simplify( (f - A*x .* x) ./ x );

f_fh = matlabFunction(f,'vars',{x,mu});             % later used for ODE solver
J_fh = matlabFunction(jacobian(f,x),'vars',{x,mu}); % compute trace and determinant of Jacobian
f1_fh = matlabFunction(f(1),'vars',[x;mu]);         % to plot vector field
f2_fh = matlabFunction(f(2),'vars',[x;mu]);         % to plot vector field

A_ = num2cell(A');
[a11,a12,a21,a22] = deal(A_{:});
b1 = double(b(1));

% %%

LimX = [-1,5];
LimY = [-1,5];
hX = 0.05;
hY = 0.05;

B = 0.3;
E = linspace(-B,B,101);
[x_,y_] = meshgrid(E,E);

mu = 1.5;
b2 = mu;

EqX = [-b1/a11 ; 0];
EqY = [0 ; -b2/a22];
EqP = -A\double(subs(b));

odefun = @(t,x) [f1_fh(x(1),x(2),mu);f2_fh(x(1),x(2),mu)];

plotpp(odefun, 'arrowSize', 12,...
    'plotQuiver', false, ...
    'axisMarginRatio', 0.5,...
    'plotNonSaddleTrajectory', false, ...
    'arrowDensity', 0.2,...
    'xlim', LimX, 'ylim', LimY,...
    'lineColor',[1,0,0])
hold on, box on, grid on;

fill([2 2 -1 -1 0 0]*LimX(2),[0 -1 -1 2 2 0]*LimY(2),[0,0,0],'FaceAlpha',0.1)

% Compute the values of the vector field in a given number of grid points.
[xx,yy] = meshgrid(LimX(1):hX:LimX(2),LimY(1):hY:LimY(2));
f1_val = f1_fh(xx,yy,mu);
f2_val = f2_fh(xx,yy,mu);

% Normalize the vectors of the vector field.
r = sqrt(f1_val.^2 + f2_val.^2);
f1_val = f1_val ./ r;
f2_val = f2_val ./ r;

streamslice(xx,yy,f1_val,f2_val);
title(str,'Interpreter','latex')