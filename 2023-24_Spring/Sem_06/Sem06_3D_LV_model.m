%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. October 25. (2023a)
%

%{

pcz_num2str_multiline(A,b,'round',20,'format','%20.17f')

%}

eta_c = 730 \ (26771 + sqrt(80211121));
nu_c = 117 - 2*eta_c;

eta = eta_c + 1/10;
nu = nu_c + 1/1000;

eta = eta_c + 1/20;
nu = nu_c + 1/2000;

eta = eta_c + 1/100;
nu = nu_c + 1/10000;

eta = eta_c + 1/1000;
nu = nu_c + 1/100000;

% eta = 0.253814;
% nu = 32 - eta;

A = -[
    23  9 eta
    37 16  nu
    9   7  16
    ];
 
b = -sum(A,2);

x = sym('x',[3,1]);
f = diag(x)*(A*x+b);
f_fh = @(t,x) diag(x)*(A*x+b);
f1_fh = matlabFunction(f(1),'vars',x);
f2_fh = matlabFunction(f(2),'vars',x);
f3_fh = matlabFunction(f(3),'vars',x);

sol = solve(f,x);
Eq3 = double([sol.x1,sol.x2,sol.x3]');
Eq3(:,any(Eq3 < 0,1)) = [];

Lims = [5,6,3];

xyz = 'xyz';
Coords = [
    1 2
    1 3 
    2 3
    ]';

term_event_2D = @(t,x) hp_ode_terminal_event_ball(t,x,[0;0],norm(Lims));
term_event_3D = @(t,x) hp_ode_terminal_event_ball(t,x,[0;0;0],norm(Lims));

initial_conditions = {
% (x,y)
  [ 0.1 6 ];
% (x,z)
  [ 5 0.1 ];
% (y,z)
  [ 0.030 3 ];
    };


fig = figure(312);
Tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
Ax1 = nexttile;
hold on, grid on, box on;


for dim = 1:size(Coords,2)
    coords = Coords(:,dim);
    A_ = A(coords,coords);
    b_ = b(coords);    
    Lims_ = Lims(coords);

    %%% 
    % Generate anonymous functions
    % 
    x = sym('x',[2,1]);
    f = diag(x)*(A_*x+b_);
    f_2D_fh = @(t,x) diag(x)*(A_*x+b_);
    J_fh = matlabFunction(jacobian(f,x),'vars',{x});
    f_c1_fh = matlabFunction(f(1),'vars',x);
    f_c2_fh = matlabFunction(f(2),'vars',x);
    
    %%%
    % Compute equilibria
    %
    Eq0 = [0;0];
    Eq1 = -A_\b_;
    EqX = [-b_(1)/A_(1,1) ; 0];
    EqY = [0 ; -b_(2)/A_(2,2)];
    Eq = [ Eq0 Eq1 EqX EqY ];
    Eq(:,any(Eq < 0,1)) = [];
        
    %%%
    %
    Resolution = 31;
    [x1,x2] = meshgrid(linspace(0,Lims(coords(1)),Resolution),linspace(0,Lims(coords(2)),Resolution));
    
    f1 = f_c1_fh(x1,x2);
    f2 = f_c2_fh(x1,x2);
    norm_f = sqrt(f1.^2 + f2.^2);
    f1 = f1 ./ norm_f;
    f2 = f2 ./ norm_f;

    O = x2*0;
    o = Eq(1,:)*0;
    switch xyz(coords)
        case 'xy'
            quiver3(x1,x2,O,f1,f2,O)
            plot3(Eq(1,:),Eq(2,:),o,'.','MarkerSize',50)
        case 'xz'
            quiver3(x1,O,x2,f1,O,f2)
            plot3(Eq(1,:),o,Eq(2,:),'.','MarkerSize',50)
        case 'yz'
            quiver3(O,x1,x2,O,f1,f2)
            plot3(o,Eq(1,:),Eq(2,:),'.','MarkerSize',50)
    end
    
    xlim([0,Lims(1)]);
    ylim([0,Lims(2)]);
    zlim([0,Lims(3)]);
    
    
    for x0 = initial_conditions{dim}'
        [~,x] = ode89(f_2D_fh,[0 10],x0,odeset('Events',term_event_2D));

        O = x(:,1)*0;
        switch xyz(coords)
            case 'xy'
                Pl = plot3(x(1,1),x(1,2),0,'.','MarkerSize',25);
                plot3(x(:,1),x(:,2),O,'Color',Pl.Color,'LineWidth',1.5);
            case 'xz'
                Pl = plot3(x(1,1),0,x(1,2),'.','MarkerSize',25);
                plot3(x(:,1),O,x(:,2),'Color',Pl.Color,'LineWidth',1.5);
            case 'yz'
                Pl = plot3(0,x(1,1),x(1,2),'.','MarkerSize',25);
                plot3(O,x(:,1),x(:,2),'Color',Pl.Color,'LineWidth',1.5);
        end
    end

end


Resolution = 11;
[X,Y,Z] = ndgrid(linspace(0,Lims(1),Resolution), ...
                 linspace(0,Lims(2),Resolution), ...
                 linspace(0,Lims(3),Resolution));

f1 = f1_fh(X,Y,Z);
f2 = f2_fh(X,Y,Z);
f3 = f3_fh(X,Y,Z);
norm_f = sqrt(f1.^2 + f2.^2 + f3.^2) + eps;
f1 = f1 ./ norm_f;
f2 = f2 ./ norm_f;
f3 = f3 ./ norm_f;
% quiver3(X,Y,Z,f1,f2,f3)

T = 100;
for x0 = [
        0.1 0.1 0.1
        0.1 0.1 0.1
        0.1 2 0.001
        1 2 1.1
        2 2 1.1
        3 3 1.1 
        Lims
        [initial_conditions{1} initial_conditions{1}(:,1)*0+0.01]
        [initial_conditions{2}(:,1) initial_conditions{2}(:,1)*0+0.01 initial_conditions{2}(:,2)]
        ]'

    [~,x] = ode89(f_fh,[0 10],x0,odeset('Events',term_event_3D));
    Pl = plot3(x(1,1),x(1,2),x(1,3),'.','MarkerSize',25);
    plot3(x(:,1),x(:,2),x(:,3),'Color',Pl.Color,'LineWidth',1.5);

end

[~,x] = ode89(f_fh,[0 2000],[0.1,0.1,0.1]',odeset('Events',term_event_3D));
x0_ = x(end,:)';
[t,LC] = ode89(f_fh,[0 T],x0_,odeset('Events',@(t,x) hp_ode_terminal_event_ball(t,x,x0_,0.1,Direction=-1),'MaxStep',0.01));
if t(end) < T
    LC = [LC ; x0_'];
end
plot3(LC(:,1),LC(:,2),LC(:,3),'Color',[0 0 0],'LineWidth',5);

[~,i] = min(LC(:,3));
r = LC(i,:);
dr = LC(i+1,:) - LC(i,:);
dr = dr / norm(dr);
LC_direction_args = {r(1),r(2),r(3),dr(1),dr(2),dr(3),[0,0,0],0};
Qv = quiver3d(LC_direction_args{:});


xlabel 'x'
ylabel 'y'
zlabel 'z'

% (x,z)
view([0,0])

% (y,z)
view([90,0])

% 3D view
view([136 26])

%% Generate specific trajectories, to plot the invariant manifold

x0 = Eq3(:,Eq3(3,:) == 0 & Eq3(1,:) == 0 & Eq3(2,:) ~= 0) + [0.0001;0;0];
[~,Tr_xy] = ode89(f_fh,[0 1],x0,odeset('Events',term_event_3D));

x0 = Eq3(:,Eq3(2,:) == 0 & Eq3(3,:) == 0 & Eq3(1,:) ~= 0) + [0;0;0.01];
[~,Tr_xz] = ode89(f_fh,[0 10],x0,odeset('Events',term_event_3D));

x0 = Eq3(:,Eq3(2,:) == 0 & Eq3(1,:) == 0 & Eq3(3,:) ~= 0) + [0;0.001;0];
[~,Tr_yz] = ode89(f_fh,[0 2],x0,odeset('Events',term_event_3D));

x0 = Eq3(:,Eq3(3,:) == 0 & Eq3(1,:) == 0 & Eq3(2,:) ~= 0) + [0.001;0;0.001];
[~,x_outer] = ode89(f_fh,[0 1000],x0,odeset('Events',term_event_3D));
% plot3(x_outer(:,1),x_outer(:,2),x_outer(:,3));

x0 = [0.93;0.94;0.98];
[~,x_inner] = ode89(f_fh,[0 5000],x0,odeset('Events',term_event_3D));
% plot3(x_inner(:,1),x_inner(:,2),x_inner(:,3));

M = 5000;

% Three equilibria on the coordinate axes. Compute the normal vector of the plane defined
% by the three equilibria.
ABC = Eq3(:,sum(Eq3 == 0) == 2);
A = ABC(:,1);
uv = diff(ABC,[],2);
n = cross(uv(:,1),uv(:,2));
n = n / norm(n);

% Collect the points
All_Points = [
    ABC'
    Tr_xy
    Tr_xz
    Tr_yz
    x_outer
    x_inner(100:end,:)
    ];

% Resample the points to have a uniform density
Idx = randsample(size(All_Points,1),1);
Points = [ All_Points(Idx,:) ; zeros(M-1,3)];
All_Points(Idx,:) = [];
r = 0.05;
for Sdx = 2:M
    Norm = vecnorm(All_Points - Points(Sdx-1,:),2,2);
    All_Points(Norm < r,:) = [];

    Idx = randsample(size(All_Points,1),1);
    Points(Sdx,:) = All_Points(Idx,:);
    All_Points(Idx,:) = [];    

    if size(All_Points,1) < 1
        break
    end
end
Points(vecnorm(Points,Inf,2) == 0,:) = [];

% Project the points on the plane defined by the three equilibria.
Points_mer = (Points-A') - (Points-A') * n * n' + A';

% Rotate the projected points to the (x,y) plane
Cov = cov(Points_mer);
[S,D] = eig(Cov);
Points_2D = (Points_mer-A')*S;

% Triangulate
shp = alphaShape(Points_2D(:,[2,3]),2);
DT = shp.alphaTriangulation;

%% Plot invariant manifold with a few trajectories

Ax2 = nexttile;
hold on, grid on, box on;

Tr = trisurf(DT,Points(:,1),Points(:,2),Points(:,3));
Tr.FaceAlpha = 0.5;
shading interp

plot3(LC(:,1),LC(:,2),LC(:,3),'Color',[0 0 0],'LineWidth',5);
Qv = quiver3d(LC_direction_args{:});

for x0 = [
        0.1 0.1 0.1
        0.1 0.1 0.1
        0.1 2 0.001
        1 2 1.1
        2 2 1.1
        3 3 1.1 
        Lims
        [initial_conditions{1} initial_conditions{1}(:,1)*0+0.01]
        [initial_conditions{2}(:,1) initial_conditions{2}(:,1)*0+0.01 initial_conditions{2}(:,2)]
        ]'

    [~,x] = ode89(f_fh,[0 10],x0,odeset('Events',term_event_3D));
    Pl = plot3(x(1,1),x(1,2),x(1,3),'.','MarkerSize',25);
    plot3(x(:,1),x(:,2),x(:,3),'Color',Pl.Color,'LineWidth',1.5);
end

view([60 25])
xlabel 'x'
ylabel 'y'
zlabel 'z'
xlim([0,Lims(1)]);
ylim([0,Lims(2)]);
zlim([0,Lims(3)]);

linkprop([Ax1,Ax2],'View')
