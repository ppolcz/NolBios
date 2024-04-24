%% Kompetitív Lotka-Volterra modell
% Két egymással versengő populáció dinamikájának modellje.
% 
% (Pl. egy kovász kultúra és egy invazívabb élesztőgomba kultúra. Állítólag 
% az élesztőgombák kiszorítják a kovászt, valószínűleg azért mert az élesztőgombák 
% gyorsabban szaporodnak.)

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];
Color_X = [0 0.4470 0.7410];
Color_Y = [0.8500 0.3250 0.0980];

%%

syms t x y mu real
F_{1} = [
    x*(1 - x/2) - x*y
    y*(mu - y) - x*y
    ];
F_{2} = [
    x*(2 - 2*x) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

f = F_{1};

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

LimX = 2.5;
LimY = 2.5;
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
    'xlim', [-0.5, LimX], 'ylim', [-0.5, LimY],...
    'lineColor',[1,0,0])
hold on, box on, grid on;

fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)

% Compute the values of the vector field in a given number of grid points.
[xx,yy] = meshgrid(-0.5:hX:LimX,-0.5:hY:LimY);
f1_val = f1_fh(xx,yy,mu);
f2_val = f2_fh(xx,yy,mu);

% Normalize the vectors of the vector field.
r = sqrt(f1_val.^2 + f2_val.^2);
f1_val = f1_val ./ r;
f2_val = f2_val ./ r;

streamslice(xx,yy,f1_val,f2_val);
