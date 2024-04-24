%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. October 25. (2023a)
%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

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
J_fh = matlabFunction(jacobian(f,x),'vars',{x});
f1_fh = matlabFunction(f(1),'vars',x);
f2_fh = matlabFunction(f(2),'vars',x);
f3_fh = matlabFunction(f(3),'vars',x);

sol = solve(f,x);
Eq = double([sol.x1,sol.x2,sol.x3]');
Eq(:,any(Eq < 0,1)) = [];

XLim = [-0.5,7];
YLim = [-0.5,7];
ZLim = [-0.5,3];
termevent = @(t,x) hp_ode_terminal_event_pp3(t,x,XLim,YLim,ZLim,Eq,0.005);
opts = odeset('Events',termevent);

fig = figure(12);
delete(fig.Children)
ax = axes(fig);
hold on, grid on, box on;

if true
    %% Vektormezo (x,y)
    
    [x,y] = meshgrid(linspace(XLim(1),XLim(2),101),linspace(YLim(1),YLim(2),101));
    
    f1 = f1_fh(x,y,y*0);
    f2 = f2_fh(x,y,y*0);
    
    r = sqrt(f1.^2 + f2.^2);
    f1 = f1 ./ r;
    f2 = f2 ./ r;
    
    [vert,avert] = streamslice(x,y,f1,f2);
    for i = 1:numel(vert)
        plot(vert{i}(:,1),vert{i}(:,2),'Color',Color_1);
    end
    for i = 1:numel(avert)
        plot(avert{i}(:,1),avert{i}(:,2),'Color',Color_1);
    end
    
    %% Vektormezo (x,z)
    
    [x,z] = meshgrid(linspace(XLim(1),XLim(2),101),linspace(ZLim(1),ZLim(2),101));
    
    f1 = f1_fh(x,x*0,z);
    f3 = f3_fh(x,x*0,z);
    
    r = sqrt(f1.^2 + f3.^2);
    f1 = f1 ./ r;
    f3 = f3 ./ r;
    
    [vert,avert] = streamslice(x,z,f1,f3);
    for i = 1:numel(vert)
        plot3(vert{i}(:,1),0*vert{i}(:,2),vert{i}(:,2),'Color',Color_2);
    end
    for i = 1:numel(avert)
        plot3(avert{i}(:,1),0*avert{i}(:,2),avert{i}(:,2),'Color',Color_2);
    end
    
    %% Vektormezo (y,z)
    
    [y,z] = meshgrid(linspace(YLim(1),YLim(2),101),linspace(ZLim(1),ZLim(2),101));
    
    f2 = f2_fh(y*0,y,z);
    f3 = f3_fh(y*0,y,z);
    
    r = sqrt(f2.^2 + f3.^2);
    f2 = f2 ./ r;
    f3 = f3 ./ r;
    
    [vert,avert] = streamslice(y,z,f2,f3);
    for i = 1:numel(vert)
        plot3(0*vert{i}(:,2),vert{i}(:,1),vert{i}(:,2),'Color',Color_5);
    end
    for i = 1:numel(avert)
        plot3(0*avert{i}(:,2),avert{i}(:,1),avert{i}(:,2),'Color',Color_5);
    end
end

%%

Pl = plot3(Eq(1,:),Eq(2,:),Eq(3,:),'.k','MarkerSize',25);

r = 0.02;
for eq = Eq

    [S,lambda] = eig(J_fh(eq),'vector');

    if any(abs(imag(lambda)) > 0)
        s1 = S(:,imag(lambda) == 0);

        % quiver3(1,1,1,s1(1),s1(2),s1(3))
        S = s1;
    end    

    Directions = [S , -S];
    for v = Directions

        x0 = eq + v*r;
        try
            [t,x] = ode45(f_fh,[0,100],x0,opts);
            plot3(x(:,1),x(:,2),x(:,3),'r','LineWidth',2)
        catch
            eq,v,x0
        end

        try
            [t,x] = ode45(f_fh,[0,-100],x0,opts);
            plot3(x(:,1),x(:,2),x(:,3),'r','LineWidth',2)
        catch
            eq,v,x0
        end

    end
end

xlim(XLim);
ylim(YLim);
zlim(ZLim);

%% Generate specific trajectories, to plot the invariant manifold

Lims = [XLim ; YLim ; ZLim];
Lims = Lims(:,2)';
term_event_3D = @(t,x) hp_ode_terminal_event_ball(t,x,[0;0;0],norm(Lims));

x0 = Eq(:,Eq(3,:) == 0 & Eq(1,:) == 0 & Eq(2,:) ~= 0) + [0.0001;0;0];
[~,Tr_xy] = ode89(f_fh,[0 1],x0,odeset('Events',term_event_3D));

x0 = Eq(:,Eq(2,:) == 0 & Eq(3,:) == 0 & Eq(1,:) ~= 0) + [0;0;0.01];
[~,Tr_xz] = ode89(f_fh,[0 10],x0,odeset('Events',term_event_3D));

x0 = Eq(:,Eq(2,:) == 0 & Eq(1,:) == 0 & Eq(3,:) ~= 0) + [0;0.001;0];
[~,Tr_yz] = ode89(f_fh,[0 2],x0,odeset('Events',term_event_3D));

x0 = Eq(:,Eq(3,:) == 0 & Eq(1,:) == 0 & Eq(2,:) ~= 0) + [0.001;0;0.001];
[~,x_outer] = ode89(f_fh,[0 1000],x0,odeset('Events',term_event_3D));
% plot3(x_outer(:,1),x_outer(:,2),x_outer(:,3));

x0 = [0.93;0.94;0.98];
[~,x_inner] = ode89(f_fh,[0 5000],x0,odeset('Events',term_event_3D));
% plot3(x_inner(:,1),x_inner(:,2),x_inner(:,3));

M = 5000;

% Three equilibria on the coordinate axes. Compute the normal vector of the plane defined
% by the three equilibria.
ABC = Eq(:,sum(Eq == 0) == 2);
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

Tr = trisurf(DT,Points(:,1),Points(:,2),Points(:,3));
Tr.FaceAlpha = 0.5;
shading interp

T = 100;
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

for x0 = [
        0.1 0.1 0.1
        0.1 0.1 0.1
        0.1 2 0.001
        1 2 1.1
        2 2 1.1
        3 3 1.1 
        Lims
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

