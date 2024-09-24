

m = 1;
l = 1;
g = 10;
b = 0.03;

% Energy function
V_fh = @(x1,x2) m*g*l*(1-cos(x1)) + 1/2*m*l^2*x2.^2;

% Key-value pairs
args.V_Levels = m*g*l * (0:0.5:5);
args.X1Lims = [-5*pi,5*pi];  % Should be an array of dimension 1x2
args.X1Res  = 300;           % Should be an array of dimension 1x1 (i.e., scalar)
args.X2Lims = [-10,10];
args.X2Res  = 300;
args.V_threshold = Inf;

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

% Prepare data for the surface plot of the energy function
[x1_mesh,x2_mesh] = meshgrid(...
    linspace(args.X1Lims(1),args.X1Lims(2),args.X1Res),...
    linspace(args.X2Lims(1),args.X2Lims(2),args.X2Res));
V_mesh = V_fh(x1_mesh,x2_mesh);
V_mesh(V_mesh > args.V_threshold) = NaN;

fig = figure(123);
fig.Position(3:4) = [ 1324 540 ];
Tl = tiledlayout(1,1);

Ax1 = nexttile; hold on, grid on, box on
surf(x1_mesh,x2_mesh,V_mesh,'facealpha',0.5)
contour(x1_mesh,x2_mesh,V_mesh,args.V_Levels,'LineWidth',1.2)

shading interp
view([-8 25])
xlim(args.X1Lims)
ylim(args.X2Lims)

xlabel('$x_1 = \vartheta$ [rad]','Interpreter','latex')
ylabel('$x_2 = \dot{\vartheta}$ [rad/s]','Interpreter','latex')
zlabel('Energy','Interpreter','latex')
Ax1.TickLabelInterpreter = "latex";
Ax1.FontSize = 16;

Ax1.XTick = (-5:5)*pi;
Ax1.XTickLabel = {'$-5\pi$','$-4\pi$','$-3\pi$','$-2\pi$','$\pi$','$0$','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$'};
Ax1.TickLabelInterpreter = "latex";

exportgraphics(fig,'/home/ppolcz/Dropbox/Peti/PhD/Oktatas/11_NDS/2024/Tematika/fig/pendulum_E.png','ContentType','image')
