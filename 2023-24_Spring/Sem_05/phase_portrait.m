function [trajectories,starts,ends,SEP_set,UEP_set] = phase_portrait(odefun,jacfun,varargin)
%% Help for plotpp:
%
% Parameter definitions:
% -- odefun: (required) this function must be a 2-dimension autonomous
%            function with respect to (t, x).
% -- tspan: time span for ode computing, default value is 40;
% -- xLim: x-axis plot range, default value is [-5,5];
% -- yLim: y-axis plot range, default value is [-5,5];
% -- xPlotNum: trajectories that start/end from horizental border,
%              default value is 4;
% -- yPlotNum: trajectories that start/end from vertical border,
%              default value is 4;
% -- lineColor: color of lines and points, default is [0.1,0.1,0.1];
% -- quiverColor: color of quiver plots, default is [0.5,0.5,0.5];
%
%
% -- plotNonSaddleTrajectory: whether or not to plot non-saddle points
%                             trajcetories, default value is true;
% -- plotEPs:    whether or not to plot equilibrium points (EPs),
%                default value is true;
% -- plotArrows: whether or not to plot arrows, default value is true;
% -- plotQuiver: whether or not to plot quiver, default value is true.
% -- arrowDensity: density arrows on each trajectory, default is 2,
%                  and the minimum should be 0.1;
% -- arrowSize: Size of arrows on the trajectory, default is 9;
% -- quiverDensity: the density of quiver, default is 20
%
% -- axisMarginRatio: axis margin for solver so that the solver would not
%                     stop when the trajectory crosses the axis limitation
%                     until it crosses the axis margin border, which is 
%                     outside the axis limitation, default is 0.1;
%
% Example 1: Attractor
%     plotpp(@(t,x)[x(2);-x(1)-1*x(2)])
% 
% Example 2: The famous van der Pol equation
%     odefun = @(t, x) [x(2); -x(1) + 0.8*(1 - x(1)^2)*x(2)];
%     plotpp(odefun,'tspan', 20,...
%                           'quivercolor', [0.6,0.6,0.6],...
%                           'linecolor', [0.3,0.3,0.3])
% 
% Example 3: Saddle node manifold of the swing equation:
%     odefun = @(t, x) [x(2); 0.5 - sin(x(1)) - 0.2*x(2)];
%     plotpp(odefun, 'plotNonSaddleTrajectory', false,...
%                                 'plotQuiver', false, 'axisMarginRatio', 0.5,...
%                                 'arrowDensity', 0.1,...
%                                 'xlim', [-2*pi, 2*pi], 'ylim', [-3, 3])
%     set(gcf, 'unit', 'centimeters', 'position', [0,0,14,8])
%     set(gca, 'ticklabelinterpreter', 'latex',...
%              'xtick', -2*pi:pi:2*pi,...
%              'xticklabels', {'$$-2\pi$$','$$-\pi$$',...
%                              '$$0$$','$$\pi$$','$$\delta\mathrm{(rad)}$$'});
%     ylabel('$$\omega\mathrm{(p.u.)}$$', 'interpreter', 'latex');
%     legend({'Trajectory'}, 'Interpreter',"latex", 'location', 'northwest')

%% Input praser
default_tspan = 40;
default_xLim = [-5, 5];
default_yLim = [-5, 5];
default_xPlotNum = 4;
default_yPlotNum = 4;
default_arrowDensity = 2;
default_arrowSize = 9;
default_axisMarginRatio = 0.1;
default_quiverDensty = 20;

default_lineColor =   [0.2,0.2,0.2];
default_quiverColor = [0.6,0.6,0.6];

default_plotNonSaddleTrajectory = true;
default_plotEPs = true;
default_plotArrows = true;
default_plotQuiver = true;

p = inputParser;
addRequired(p, 'odefun');
addParameter(p, 'EP', []);
addParameter(p, 'tspan', default_tspan);
addParameter(p, 'xLim', default_xLim);
addParameter(p, 'yLim', default_yLim);
addParameter(p, 'xPlotNum', default_xPlotNum);
addParameter(p, 'yPlotNum', default_yPlotNum);
addParameter(p, 'arrowDensity', default_arrowDensity);
addParameter(p, 'arrowSize', default_arrowSize);
addParameter(p, 'axisMarginRatio', default_axisMarginRatio);
addParameter(p, 'quiverDensity', default_quiverDensty);


addParameter(p, 'plotNonSaddleTrajectory',...
                 default_plotNonSaddleTrajectory);
addParameter(p, 'plotEPs', default_plotEPs);
addParameter(p, 'plotArrows', default_plotArrows);
addParameter(p, 'plotQuiver', default_plotQuiver);

parse(p, odefun, varargin{:});
r = p.Results;

tspan = r.tspan;
xLim = r.xLim;
yLim = r.yLim;
xPlotNum = r.xPlotNum;
yPlotNum = r.yPlotNum;
arrowDensity = r.arrowDensity;
    
%% Initialization
% ax = axes(figure);
% set(ax, 'box', 'on','nextplot', 'add',...
%     'xgrid', 'on', 'ygrid', 'on',...
%     'xlim', xLim, 'ylim', yLim);

x0_set_forward = [];
x0_set_backward = [];

EP_disturb = abs(1E-3*(xLim(1) - xLim(2)));

% tolerance border for ode
margin_ratio = r.axisMarginRatio;
x_margin = (xLim(2)-xLim(1))*margin_ratio;
y_margin = (yLim(2)-yLim(1))*margin_ratio;
xBorder = [xLim(1)-x_margin, xLim(2)+x_margin];
yBorder = [yLim(1)-y_margin, yLim(2)+y_margin];

f = @(x) odefun(0, x);

%%  Step 1: Plot quiver

%%  Step 2: Searching equilibrium points

% prepare for searching
Nmeshx = 5;
Nmeshy = 5;
tol = 1E-2;
x10 = linspace(xLim(1), xLim(2), Nmeshx);
x20 = linspace(yLim(1), yLim(2), Nmeshy);

% settins of solver
opt = optimoptions('fsolve','Algorithm', 'levenberg-marquardt',...
                    'Display', 'off',...
                    'FunctionTolerance', 1e-6);

% prepare for searching
if ~isempty(r.EP)
    EP_set = r.EP;
    EP_count = size(r.EP,2);


    ldx = r.xLim(1) < EP_set(1,:) & EP_set(1,:) < r.xLim(2) & r.yLim(1) < EP_set(2,:) & EP_set(2,:) < r.yLim(2);
    EP_set = EP_set(:,ldx);
    EP_count = sum(ldx);
else
    EP_set = [];
    EP_count = 0;

    % begin searching
    for cnt_x10 = 1:Nmeshx
        for cnt_x20 = 1:Nmeshy
            % search the EP from each meshgrid point
            pending_EP = fsolve(f, [x10(cnt_x10); x20(cnt_x20)], opt);
            if (pending_EP(1)>xLim(1)) && (pending_EP(1)<xLim(2)) &&...
                (pending_EP(2)>yLim(1)) && (pending_EP(2)<yLim(2))
                % make sure the EPs are inside the plot region
                addIntoSet = true;
                % check if there is the same EP
                if ~isempty(EP_set)
                    EP_err = EP_set-pending_EP;
                    % this is realized by checking the distance
                    if any(vecnorm(EP_err,2,2)<tol)
                       addIntoSet = false;
                    end
                end
            else
                addIntoSet = false;
            end
            if addIntoSet
                EP_set = [EP_set , pending_EP];
                EP_count = EP_count + 1;
            end
        end
    end
end

% send EP_set to the base workspace
assignin('base', 'EP_set', EP_set);

%% Step 3: Checking the type of equilibrium points
if EP_count == 0
else
    SEP_set = []; UEP_set = [];
    vec_len = EP_disturb;
    for cnt_EP = 1:EP_count
        %% calc the Jacobian
        EP = EP_set(:,cnt_EP);
        jcb_mat = jacfun(EP);
        [V, D] = eig(jcb_mat, 'vector');

        D1_real = real(D(1));
        D2_real = real(D(2));
        D_imag = imag(D);
        
        %% decide the EP type from Jacobian
        if (D1_real<0) && (D2_real<0)
            % stable equilibrium point
            SEP_set = [SEP_set , EP];
            if (r.plotNonSaddleTrajectory)
            if (D_imag) == 0
                % stable node (need four trajectories)
                x0_set_backward = [x0_set_backward , EP + V(:,1)*vec_len];
                x0_set_backward = [x0_set_backward , EP - V(:,1)*vec_len];
                x0_set_backward = [x0_set_backward , EP + V(:,2)*vec_len];
                x0_set_backward = [x0_set_backward , EP - V(:,2)*vec_len];
            else 
                % stable focus (only need one trajecotry)
                x0_set_backward = [x0_set_backward , EP + [1,0]*vec_len];
                % x0_set_backward = [x0_set_backward , EP + [0,1]*vec_len];
                % x0_set_backward = [x0_set_backward , EP - [1,0]*vec_len];
                % x0_set_backward = [x0_set_backward , EP - [0,1]*vec_len];
            end
            end
        elseif (D1_real>0) && (D2_real>0)
            UEP_set = [UEP_set , EP];
            if (r.plotNonSaddleTrajectory)
            if (D_imag) == 0
                % unstable node (need four trajectories)
                x0_set_forward = [x0_set_forward , EP + V(:,1)*vec_len];
                x0_set_forward = [x0_set_forward , EP - V(:,1)*vec_len];
                x0_set_forward = [x0_set_forward , EP + V(:,2)*vec_len];
                x0_set_forward = [x0_set_forward , EP - V(:,2)*vec_len];
            else 
                % unstable focus (only need one trajecotry)
                x0_set_forward = [x0_set_forward , EP + [1,0]*vec_len];
            end
            end
        elseif (D1_real == 0) || (D2_real == 0)
            UEP_set = [UEP_set , EP];
            % center (only need one trajecotry)
            x0_set_forward = [x0_set_forward , EP + [1,0]*vec_len];
        else
            UEP_set = [UEP_set , EP];
            % saddle node (need four trajectories)
            if (D1_real>0)
                x0_set_forward = [x0_set_forward , EP + V(:,1)*vec_len];
                x0_set_forward = [x0_set_forward , EP - V(:,1)*vec_len];
                x0_set_backward = [x0_set_backward , EP + V(:,2)*vec_len];
                x0_set_backward = [x0_set_backward , EP - V(:,2)*vec_len];
            else
                x0_set_backward = [x0_set_backward , EP + V(:,1)*vec_len];
                x0_set_backward = [x0_set_backward , EP - V(:,1)*vec_len];
                x0_set_forward = [x0_set_forward , EP + V(:,2)*vec_len];
                x0_set_forward = [x0_set_forward , EP - V(:,2)*vec_len];
            end
        end
    end
end

%%  Step 5: Solve ODEs and plot phase portrait
% control parameters
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-7,...
              'Events', @(t,x) myevent(t,x,xBorder,yBorder));
starts = []; ends = [];

trajectories = {};

forward_plotNum = size(x0_set_forward,2);
for cnt = 1:forward_plotNum
    [~,xsol] = ode45(odefun, [0 tspan], x0_set_forward(:,cnt), opts);
    trajectories{end+1} = xsol;
    if (r.plotArrows)
    [starts, ends] = calcSpecGrad(starts, ends,...
                                  xsol(:,1), xsol(:,2),...
                                  xLim, yLim, arrowDensity);
    end
end

backward_plotNum = size(x0_set_backward,2);
for cnt = 1:backward_plotNum
    [~,xsol] = ode45(odefun, [0 -tspan], x0_set_backward(:,cnt), opts);
    trajectories{end+1} = xsol;
    % plot(ax, xsol(:,1), xsol(:,2),...
    %          'color', cline, 'linewidth', 1.0);
    if (r.plotArrows)
    [starts, ends] = calcSpecGrad(starts, ends,...
                                  flip(xsol(:,1)), flip(xsol(:,2)),...
                                  xLim, yLim, arrowDensity);
    end
end

end

%% trigger condition for ode
function [value,isterminal,direction] = myevent(~,x,xBorder,yBorder)
    
    % if system state pass the border
    if (x(1)-xBorder(1)>0)&&(xBorder(2)-x(1)>0)&&...
       (x(2)-yBorder(1)>0)&&(yBorder(2)-x(2)>0)
        value = 1;
    else
        value = -1;
    end
    
    isterminal=1;
    direction= -1;
end

%% functions from phasePortraitPlot
function [starts, ends] = calcSpecGrad(starts, ends, x, y,...
                                       xLim, yLim, arrowDensity)
    xscale = xLim(2) - xLim(1);
    yscale = yLim(2) - yLim(1);
    ds = sqrt(gradient(x/xscale).^2 + gradient(y/yscale).^2);
    s = cumsum(ds);
    num = ceil(s(end)*arrowDensity);
    sTick = s(end)/num*((1:num)-0.5);
    
    for cnt = 1:numel(sTick)
        pos = find(s >= sTick(cnt), 1);
        loc = [x(pos), y(pos)];
        if (loc(1)>xLim(1))&&(loc(1)<xLim(2)&&...
            loc(2)>yLim(1))&&(loc(2)<yLim(2)) 
            if isempty(starts)
                starts = [starts; loc];
                ends = [ends; [x(pos+1), y(pos+1)]];
            elseif all(vecnorm((loc-starts)./[xscale,yscale],2,2)>0.1)
                starts = [starts; loc];
                ends = [ends; [x(pos+1), y(pos+1)]];
            end
        end
    end
end

%% function for plotting arrows
function arrows(starts, ends, arrowSizePix, varargin)
    [hei, ~] = size(starts);
    for cnt = 1:hei
        arrow(starts(cnt,:), ends(cnt,:), arrowSizePix, 'BaseAngle', 60,...
             varargin{:});
    end
end