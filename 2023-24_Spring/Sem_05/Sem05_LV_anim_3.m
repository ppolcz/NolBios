%% Kompetitív Lotka-Volterra modell
% Két egymással versengő populáció dinamikájának modellje.
% 
% (Pl. egy kovász kultúra és egy invazívabb élesztőgomba kultúra. Állítólag 
% az élesztőgombák kiszorítják a kovászt, valószínűleg azért mert az élesztőgombák 
% gyorsabban szaporodnak.)

%%

xe = 1;
ye = 2;

syms t x y mu eta K real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
assumeAlso(K >= 1)
assumeAlso(K <= 5)
f = [
    x*( ye*(1-x/K) - y*(1-xe/K) )
    y*( xe*(1-y/K) - x*(1-ye/K) )
    ];
x = [x;y];

sol = solve(f);
Eq = [
    sol.x sol.y
    ]';

f_fh = matlabFunction(f,'vars',{x,K});             % later used for ODE solver
J_fh = matlabFunction(jacobian(f,x),'vars',{x,K}); % compute trace and determinant of Jacobian
f1_fh = matlabFunction(f(1),'vars',[x;K]);         % to plot vector field
f2_fh = matlabFunction(f(2),'vars',[x;K]);         % to plot vector field

% A_ = num2cell(A');
% [a11,a12,a21,a22] = deal(A_{:});
% b1 = double(b(1));

str = "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";

%%

LimX = 5;
LimY = 5;
hX = 0.05;
hY = 0.05;

B = 0.3;
E = linspace(-B,B,101);
[x_,y_] = meshgrid(E,E);

fig = figure(1);
fig.Position(3:4) = [1500 450];

mu_val = 1.3;

ind = 0;
for K_val = 1.5:0.1:4.5
    ind = ind + 1;
        
    odefun = @(t,x) f_fh(x,K_val);
    jacfun = @(x) J_fh(x,K_val);
    
    [trajectories,starts,ends,SEP_set,UEP_set] = phase_portrait(odefun,jacfun, ...
        'EP',double(subs(Eq,[K],[K_val])), ...
        'xlim', [-0.5, LimX],'ylim', [-0.5, LimY], ...
        'plotNonSaddleTrajectory',true, ...
        'tspan',500);
    
    % Compute the values of the vector field in a given number of grid points.
    [xx,yy] = meshgrid(-0.5:hX:LimX,-0.5:hY:LimY);
    f1_val = f1_fh(xx,yy,K_val);
    f2_val = f2_fh(xx,yy,K_val);
    
    % Normalize the vectors of the vector field.
    r = sqrt(f1_val.^2 + f2_val.^2);
    f1_val = f1_val ./ r;
    f2_val = f2_val ./ r;
        
    Tl = tiledlayout(1,3);
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
    plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    xlim([-0.5,LimX])
    ylim([-0.5,LimY])
    
    fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    
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
    plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    xlim([-0.5,LimX])
    ylim([-0.5,LimY])
    title(str + ", $K = " + num2str(K_val) + "$",'Interpreter','latex')
    
    fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
        
    streamslice(xx,yy,f1_val,f2_val);

    nexttile, hold on, box on, grid on
    streamslice(xx,yy,f1_val,f2_val);
    if ~isempty(SEP_set)
        plot(SEP_set(1,:),SEP_set(2,:),'r.','MarkerSize',25);
    end
    plot(UEP_set(1,:),UEP_set(2,:),'ro','MarkerSize',7,'LineWidth',2);
    fill([2 2 -1 -1 0 0]*LimX,[0 -1 -1 2 2 0]*LimY,[0,0,0],'FaceAlpha',0.1)
    xlim([-0.5,LimX])
    ylim([-0.5,LimY])

    fimplicit(@(x,y) f1_fh(x,y,K_val),[-0.5,LimX,-0.5,LimY]);
    fimplicit(@(x,y) f2_fh(x,y,K_val),[-0.5,LimX,-0.5,LimY]);


    drawnow
    % exportgraphics(fig,sprintf('Results/%03d.jpg',ind))
end
