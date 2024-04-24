%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. October 25. (2023a)
%

r = 1.5;
a = 5;
b = 18;
D = 2;
c = 4;

K_Hopf = b*(a*D+c)/(a*D-c);

syms t x y K real
f = [
    r*x*(1 - x/K) - y*(a*x)/(b+x)
    -c*y + D*y*(a*x)/(b+x)
    ];
x = [x;y];

f_fh = matlabFunction(f,'vars',{x,K});             % later used for ODE solver
J_fh = matlabFunction(jacobian(f,x),'vars',{x,K}); % compute trace and determinant of Jacobian
f1_fh = matlabFunction(f(1),'vars',[x;K]);         % to plot vector field
f2_fh = matlabFunction(f(2),'vars',[x;K]);         % to plot vector field

LimX = K_Hopf*1.5;
LimY = K_Hopf;

% This is not known ...
Eq1_fh = @(K) [
    b*c/(a*D-c)
    b*D*r*(-K*c-b*c+K*a*D)/K/(a*D-c)^2
    ];

K_vals = unique([...
    linspace(K_Hopf/2,K_Hopf*1.4,51),...         coarse grid
    linspace(K_Hopf*0.99,K_Hopf*1.01,51),...     fin grid around critical point
    ]);

% Initial guess for the first equilibrium point
Eq1_guess = Eq1_fh(K_vals(1)) + randn(2,1)*0.1;

% Equilibrium points will be stored here.
SS = [0;0]*K_vals*NaN;

% Limit cycles will be stored here.
LCs = cell(size(K_vals));

fig = figure(312);

% Loop through the possible K values
Img = imread('trace-det-ref-1.png');
for i = 1:numel(K_vals)
    K_val = K_vals(i);

    SS(:,i) = Eq1_fh(K_val);

    Tl = tiledlayout(2,2,"Padding","compact");
    
    ax = nexttile;
    hold on, box on, grid on;

    % Visualize the null-clines using function FIMPLICIT.
    fimplicit(@(x,y) f1_fh(x,y,K_val),[0 LimX],'Color',COL.Color_2,'LineWidth',1.5);  % <---------- FIMPLICIT
    fimplicit(@(x,y) f2_fh(x,y,K_val),[0 LimY],'Color',COL.Color_2,'LineWidth',1.5);
    plot(SS(1,i),SS(2,i),'.k','MarkerSize',15);
    
    % Compute the values of the vector field in a given number of grid points.
    [xx,yy] = meshgrid(linspace(0,LimX,50),linspace(0,LimY,50));
    f1_val = f1_fh(xx,yy,K_val);
    f2_val = f2_fh(xx,yy,K_val);
    
    % Normalize the vectors of the vector field.
    r = sqrt(f1_val.^2 + f2_val.^2);
    f1_val = f1_val ./ r;
    f2_val = f2_val ./ r;
    
    % Visualize the normalized vector field.
    % Qv = quiver(xx,yy,f1_val,f2_val,'Color',[1,1,1]*0.7);
    Ln = streamslice(xx,yy,f1_val,f2_val);
    for l = Ln'
        l.Color = COL.Color_Gray;
    end
    
    % Solve the ODE from a given initial condition and visualize its trajectory.
    [~,xx] = ode89(@(t,x) f_fh(x,K_val),[0,100],[1;1]);
    plot(xx(:,1),xx(:,2),'Color',COL.Color_1,'LineWidth',1.5);

    if K_val >= K_Hopf
        % If K is above the critical value, compute and store the limit cycle.

        term_event = @(t,x) hp_ode_terminal_event_ball(t,x,xx(end,:)',0.5,"Direction",-1);
        opts = odeset('Events',term_event,'MaxStep',0.01);
        [~,LC] = ode89(@(t,x) f_fh(x,K_val),[0,100],xx(end,:)',opts);
        LC = [ LC ; LC(1,:) ];

        plot(LC(:,1),LC(:,2),'r','LineWidth',2);
        LCs{i} = LC;
    end
    
    xlim([0 LimX]);
    ylim([0 LimY]);

    xlabel('prey (x)')
    ylabel('predator (y)')
    title(ax,sprintf('K = %g',K_val))

    % Trace-determinant map
    ax = nexttile;
    hold on, grid on, box on;

    Tr = linspace(-3,3,101);
    Pl = plot(Tr,Tr.^2/4,'Color',COL.Color_2);
    J = J_fh(SS(:,i),K_val);
    d = det(J);
    t = trace(J);
    plot(t,d,'.k','MarkerSize',25)
    title(ax,sprintf('K = %g',K_val))

    XLim = ax.XLim;
    xgrid = linspace(XLim(1),XLim(2),10000);
    YLim = [-1,3];
    plot(XLim,[0 0],'Color',Pl.Color);
    plot([0 0],YLim,'Color',Pl.Color);
    plot(xgrid,xgrid.^2/4,'Color',Pl.Color);
    ax.XLim = XLim;
    ax.YLim = YLim;
    xlabel('trace(A)')
    ylabel('det(A)')
    
    % Equilibrium point as a function of K.
    ax = nexttile;
    hold on, grid on, box on;
    plot3(K_vals,SS(1,:),SS(2,:),'k')
    for j = 1:i
        if ~isempty(LCs{j})
            plot3(LCs{j}(:,1)*0+K_vals(j),LCs{j}(:,1),LCs{j}(:,2),'r','LineWidth',1);
        end
    end
            
    xlim(K_vals([1,end]));
    ylim([0 LimX]);
    zlim([0 LimY]);
    view([46.9156 20.3111])
    xlabel('K')
    ylabel('x')
    zlabel('y')
    
    ax = nexttile;
    imshow(Img)

    drawnow

    fprintf('%d/%d\n',i,numel(K_vals));

end


