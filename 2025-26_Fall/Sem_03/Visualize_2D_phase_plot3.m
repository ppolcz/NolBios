%%
% Author: Peter Polcz (ppolcz@gmail.com) 
% Created on 2024. September 27. (2023a)
%
function Visualize_2D_phase_plot3(f,x,p,p_vals,opts)
arguments
    f,x
    p = sym('DUMMY__');
    p_vals = 0;
    opts.PositiveSystem = false;
    opts.EqPoints = [];
    opts.XRes = 0.05;
    opts.YRes = 0.05;
    opts.XLim = [-0.5,2.5];
    opts.YLim = [-0.5,2.5];
    opts.TimeSpan = 500;
    opts.FigNr = 141;
    opts.FigDim = [1200 370];
    opts.LaTeX = @(p) "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";
    opts.RunAfter = @(Fig,ind,p_val) [];
end

sol = solve(f,x);
Eq = [
    sol.x sol.y
    ]';

if ~isempty(p)
    f_fh = matlabFunction(f,'vars',{x,p});              % later used for ODE solver
    J_fh = matlabFunction(jacobian(f,x),'vars',{x,p});  % compute trace and determinant of Jacobian
    f1_fh = matlabFunction(f(1),'vars',[x;p]);          % to plot vector field
    f2_fh = matlabFunction(f(2),'vars',[x;p]);          % to plot vector field
end

XLim = opts.XLim;
YLim = opts.YLim;
limX = XLim(1);
LimX = XLim(2);
limY = YLim(1);
LimY = YLim(2);

Fig = figure(opts.FigNr);
Fig.Position(3:4) = opts.FigDim;
ind = 0;
for p_val = p_vals
    ind = ind + 1;
    
    p_val_cell = num2cell(p_val);

    odefun = @(t,x) f_fh(x,p_val);
    jacfun = @(x) J_fh(x,p_val);

    [trajectories,starts,ends,SEP_set,UEP_set] = phase_portrait(odefun,jacfun, ...
        'EP',double(subs(Eq,p,p_val)), ...
        'xlim',XLim,'ylim',YLim, ...
        'plotNonSaddleTrajectory',true, ...
        'tspan',opts.TimeSpan);
    
    % Compute the values of the vector field in a given number of grid points.
    [xx,yy] = meshgrid(limX:opts.XRes:LimX,limY:opts.YRes:LimY);
    f1_val = f1_fh(xx,yy,p_val_cell{:});
    f2_val = f2_fh(xx,yy,p_val_cell{:});
    
    % Normalize the vectors of the vector field.
    r = sqrt(f1_val.^2 + f2_val.^2);
    f1_val = f1_val ./ r;
    f2_val = f2_val ./ r;
        
    Tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
    nexttile, hold on, box on, grid on
    for i = 1:numel(trajectories)
        x = trajectories{i};
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);
    end
    if ~isempty(starts)
        arrow(starts,ends,12,'Color','r')
    end
    if ~isempty(SEP_set)
        plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
    end
    if ~isempty(UEP_set)
        plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    end
    xlim([limX,LimX])
    ylim([limY,LimY])
    
    if opts.PositiveSystem
        fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    end
    
    nexttile, hold on, box on, grid on
    for i = 1:numel(trajectories)
        x = trajectories{i};
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);
    end
    if ~isempty(starts)
        arrow(starts,ends,12,'Color','r')
    end
    if ~isempty(SEP_set)
        plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
    end
    if ~isempty(UEP_set)
        plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    end
    xlim([limX,LimX])
    ylim([limY,LimY])
    title(opts.LaTeX(p_val),'Interpreter','latex')
    
    if opts.PositiveSystem
        fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    end

    streamslice(xx,yy,f1_val,f2_val);
    
    nexttile, hold on, box on, grid on
    streamslice(xx,yy,f1_val,f2_val);
    if ~isempty(SEP_set)
        plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
    end
    if ~isempty(UEP_set)
        plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    end
    if opts.PositiveSystem
        fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    end
    xlim(XLim)
    ylim(YLim)
    
    fimplicit(@(x,y) f1_fh(x,y,p_val_cell{:}),[XLim,YLim]);
    fimplicit(@(x,y) f2_fh(x,y,p_val_cell{:}),[XLim,YLim]);
    
    drawnow
    opts.RunAfter(Fig,ind,p_val);
end

end
