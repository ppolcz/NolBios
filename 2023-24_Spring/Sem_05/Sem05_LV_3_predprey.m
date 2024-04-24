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
f = [
    x*(-1 + 2*x) - x*y
    y*(-4 - 3*y) + 7*x*y
    ];
x = [x;y];

sol = solve(f);
Eq = double([
    sol.x sol.y
    ]');

A = double(jacobian( simplify(f ./ x) , x ));
b = simplify( (f - A*x .* x) ./ x );

f_fh = matlabFunction(f,'vars',{x,mu});             % later used for ODE solver
J_fh = matlabFunction(jacobian(f,x),'vars',{x,mu}); % compute trace and determinant of Jacobian
f1_fh = matlabFunction(f(1),'vars',[x;mu]);         % to plot vector field
f2_fh = matlabFunction(f(2),'vars',[x;mu]);         % to plot vector field

% %%

LimX = [-0.2,2];
LimY = [-0.2,2];
hX = 0.05;
hY = 0.05;

B = 0.3;
E = linspace(-B,B,101);
[x_,y_] = meshgrid(E,E);

mu = 1.5;
b2 = mu;

odefun = @(t,x) f_fh(x,mu);
jacfun = @(x) J_fh(x,mu);

[trajectories,starts,ends,SEP_set,UEP_set] = phase_portrait(odefun,jacfun, ...
    'EP',Eq, ...
    'xlim',LimX,'ylim',LimY, ...
    'tspan',40, ...
    'plotNonSaddleTrajectory', true);

fig = figure(1);
Tl = tiledlayout(1,1);
nexttile, hold on, box on, grid on
for i = 1:numel(trajectories)
    x = trajectories{i};
    plot(x(:,1),x(:,2),'r');
end
plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);

[hei,~] = size(starts);
for cnt = 1:hei
    arrow(starts(cnt,:),ends(cnt,:),8,'BaseAngle',60,'Color','red');
end

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

% fimplicit(@(x,y) f1_fh(x,y,mu),[-0.5,LimX,-0.5,LimY]);
% fimplicit(@(x,y) f2_fh(x,y,mu),[-0.5,LimX,-0.5,LimY]);

xlim(LimX)
ylim(LimY)
