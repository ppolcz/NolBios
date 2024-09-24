function [ret] = hp_animate_pendulum(t_sol,x_sol,V_fh,u_fh,args)
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. September 28. (2023a)
%
% Possible signatures:
% hp_animate_pendulum(t_sol,x_sol,V_fh,u_fh,...
%     'V_Levels',0:0.1:5,...
%     'X1Lims',[-6,6])
arguments

    % Positional arguments
    t_sol (:,1)
    x_sol (:,2)
    V_fh
    u_fh = @(t) t*0;

    % Key-value pairs
    args.V_Levels = 0:0.5:5
    args.X1Lims (1,2) = [-5*pi,5*pi]  % Should be an array of dimension 1x2
    args.X1Res  (1,1) = 300           % Should be an array of dimension 1x1 (i.e., scalar)
    args.X2Lims (1,2) = [-10,10]
    args.X2Res  (1,1) = 300
    args.V_threshold = Inf

    args.TeXInput = '0'

    args.Record = false;
    args.ExportName = 'Pendulum'
end

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

% Compute energy along the trajectory
V_xx = V_fh(x_sol(:,1), x_sol(:,2));

T = t_sol(end);
fs = 20;
Ts = 1/fs;

t_lin = linspace(0,T,T*fs);
x_lin = interp1(t_sol,x_sol,t_lin);
V_lin = V_fh(x_lin(:,1),x_lin(:,2));

phi = x_lin(:,1);
u = u_fh(t_lin);

fig = figure(123);
fig.Position(3:4) = [ 1324 740 ];
Tl = tiledlayout(2,4);

Res = 100;
Ax3 = polaraxes(Tl); hold on, grid on, box on
IngaPl = polarplot([phi(1),phi(1)],[0,0.9],'LineWidth',10);
UPl = polarplot(linspace(0,u(1),Res)+pi,zeros(1,Res)+0.5,'LineWidth',10);
polarplot(0,0,'k.','MarkerSize',40)
title({'Illustrating the movement of the pendulum',['$ml\ddot{y} + bl\dot{y} + mg\sin{y}=' args.TeXInput '$']},'Interpreter','latex')

Ax3.ThetaZeroLocation = 'bottom';
Ax3.RLim = [0,1];
Ax3.RTick = [];
Ax3.RTickLabel = {};

Ax1 = nexttile([1,2]); hold on, grid on, box on
surf(x1_mesh,x2_mesh,V_mesh,'facealpha',0.5)
plot(x_lin(1,1), x_lin(1,2),'.k','markersize',15)
plot3(x_lin(1,1), x_lin(1,2),V_xx(1),'.k','markersize',15)
Pl2 = plot(x_lin(:,1), x_lin(:,2),'Color',Color_2,'LineWidth',2);
Pl3 = plot3(x_lin(:,1), x_lin(:,2),V_lin,'Color',Color_3,'LineWidth',2);
contour(x1_mesh,x2_mesh,V_mesh,args.V_Levels,'LineWidth',1.2)

Ax1.XTick = (-5:5)*pi;
Ax1.XTickLabel = {'$-5\pi$','$-4\pi$','$-3\pi$','$-2\pi$','$\pi$','$0$','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$'};
Ax1.TickLabelInterpreter = "latex";

shading interp
view([-8 25])
xlim(args.X1Lims)
ylim(args.X2Lims)

xlabel('$x_1 = y$ [rad]','Interpreter','latex')
ylabel('$x_2 = \dot{y}$ [rad/s]','Interpreter','latex')
Ax1.FontSize = 16;

Ax2 = nexttile; hold on, grid on, box on
plot(t_sol,V_xx), grid on
title('$E(t) := E(x(t))$', 'interpreter', 'latex','fontsize',22)
xlabel('time $t$ [s]', 'interpreter', 'latex','fontsize',22)
Ax2.TickLabelInterpreter = "latex";
Ax2.FontSize = 16;
YLim = Ax2.YLim;
Pl_Vert = plot([t_lin(1),t_lin(1)],YLim,'Color',Color_2,'LineWidth',1.5);
ylim(YLim)

AxBottom = nexttile([1,4]); hold on; grid on; box on
contour(x1_mesh,x2_mesh,V_mesh,args.V_Levels,'LineWidth',1.2)
PlO = plot(x_sol(:,1),x_sol(:,2),'Color',Color_2,'LineWidth',2);
plot(x_sol(1,1),x_sol(1,2),'Color',Color_2,'MarkerSize',25);
Font = {'Interpreter','latex','FontSize',16};
xlabel('$\theta$ [rad]',Font{:});
ylabel('$\dot\theta$ [rad/s]',Font{:});

AxBottom.XTick = (-5:5)*pi;
AxBottom.XTickLabel = {'$-5\pi$','$-4\pi$','$-3\pi$','$-2\pi$','$\pi$','$0$','$\pi$','$2\pi$','$3\pi$','$4\pi$','$5\pi$'};
AxBottom.TickLabelInterpreter = "latex";
AxBottom.FontSize = 16;


%% 

drawnow
pause(0.5)
exportgraphics(fig,string(args.ExportName) + ".png","ContentType","image")
pause(0.5)

tic;
for i = 2:numel(phi)
    Pl2.XData = x_lin(1:i,1);
    Pl2.YData = x_lin(1:i,2);
    Pl3.XData = x_lin(1:i,1);
    Pl3.YData = x_lin(1:i,2);
    Pl3.ZData = V_lin(1:i);

    Pl_Vert.XData = [t_lin(i),t_lin(i)];

    IngaPl.ThetaData = [phi(i),phi(i)];
    UPl.ThetaData = linspace(0,u(i),Res)+pi;
    if u(i) < 0
        UPl.Color = Color_2;
    else
        UPl.Color = Color_5;
    end

    PlO.XData = x_lin(1:i,1);
    PlO.YData = x_lin(1:i,2);

    drawnow
    if args.Record
        Video_Frames(i) = getframe(fig); 
    end

    % if mod(i,25) == 0
    %     keyboard
    % end

    dt = toc;
    tic;
    pause(Ts - dt);
end

if args.Record
    Video_Frames(1) = [];
    video = VideoWriter(string(args.ExportName) + ".avi", 'Uncompressed AVI');
    video.FrameRate = fs;
    open(video)
    writeVideo(video, Video_Frames);
    close(video);
end

end
