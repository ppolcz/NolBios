%%
% Author: Peter Polcz (ppolcz@gmail.com) 
% Created on 2024. September 27. (2023a)
%
function Visualize_2D_phase_plot2(f,x,p,p_vals,opts)
arguments
    f,x
    p = sym('DUMMY__');
    p_vals = 0;
    opts.EqPoints = [];
    opts.XRes = 0.05;
    opts.YRes = 0.05;
    opts.XLim = [-0.5,2.5];
    opts.YLim = [-0.5,2.5];
    opts.TimeSpan = 500;
    opts.FigNr = 141;
    opts.FigDim = [830 370];
    opts.LaTeX = @(p) "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";
    opts.RunAfter = @(Fig,ind,p_val) [];
    opts.PlotDirections = 1;
    opts.ShadeOrthants = 0;
end

sol = solve(f,x);
Eq = [
    sol.x1 sol.x2
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

    [trajectories,starts,ends,SEP_set,UEP_set,limit_cycles] = phase_portrait(odefun,jacfun, ...
        'EP',double(subs(Eq,p,p_val)), ...
        'xlim',XLim,'ylim',YLim, ...
        'plotNonSaddleTrajectory',true, ...
        'tspan',opts.TimeSpan);

    %% Generate title and directions
    % in which we present the eigenvalues of the equilibria

    EP = [SEP_set , UEP_set];
    [~,Idx] = sort(min(abs(EP),[],1) + eps*vecnorm(EP,2) + 1e3*eps*abs(EP(2,:)));
    EP = EP(:,Idx);

    e = zeros(2,0);
    s = zeros(4,0);

    u = 1;
    txt = cell(1,size(EP,2));
    for xs = EP
        [S,D] = eig(jacfun(xs));
        lambda = diag(D);
        if abs(imag(lambda(1))) > 0
            txt{u} = sprintf('Eq%d: $%.2f %s %.2fj$',u,real(lambda(1)),'\pm',imag(lambda(1)));
        else
            e = [e xs xs xs xs];
            s = [s S -S];
            txt{u} = sprintf('Eq%d: $%.2f, %.2f$', u,sort(lambda));
        end
        u = u + 1;
    end
    
    txt2 = cell(1,ceil(numel(txt)/2));
    txt2{end} = txt{end};
    for u = 1:2:numel(txt)-1
        txt2{(u+1)/2} = strjoin(txt(u:u+1),', ');
    end

    %%

    % Compute the values of the vector field in a given number of grid points.
    [xx,yy] = meshgrid(limX:opts.XRes:LimX,limY:opts.YRes:LimY);
    f1_val = f1_fh(xx,yy,p_val_cell{:});
    f2_val = f2_fh(xx,yy,p_val_cell{:});
    
    % Normalize the vectors of the vector field.
    r = sqrt(f1_val.^2 + f2_val.^2);
    f1_val = f1_val ./ r;
    f2_val = f2_val ./ r;
        
    %% Plot 1

    Tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
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
    
    if opts.ShadeOrthants
        fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    end
        
    streamslice(xx,yy,f1_val,f2_val);
    
    %% Plot 2

    nexttile, hold on, box on, grid on
    streamslice(xx,yy,f1_val,f2_val);
    if ~isempty(SEP_set)
        plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
    end
    if ~isempty(UEP_set)
        plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    end
    if opts.ShadeOrthants
        fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    end
    xlim(XLim)
    ylim(YLim)
    
    fimplicit(@(x,y) f1_fh(x,y,p_val_cell{:}),[XLim,YLim]);
    fimplicit(@(x,y) f2_fh(x,y,p_val_cell{:}),[XLim,YLim]);

    if opts.PlotDirections > 0
        ss = s*opts.PlotDirections;
        q = quiver(e(1,:),e(2,:),ss(1,:),ss(2,:));
        q.ShowArrowHead = "off";
        q.LineWidth = 2;
        q.AutoScale = "off";
        q.Color = 'red';
    end

    for i = 1:numel(limit_cycles)
        x = limit_cycles{i};
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);
    end


    title(txt2,'Interpreter','latex')

    drawnow
    opts.RunAfter(Fig,ind,p_val);
end

end
