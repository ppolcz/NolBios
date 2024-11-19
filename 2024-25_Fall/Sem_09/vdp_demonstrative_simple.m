%%
%  Author: Peter Polcz (ppolcz@gmail.com)
%  Created on November 19. (2023a)
%

%% The Van der Pol system
% 
% We consider the time-reversed Van der Pol oscillator system. The original system has a
% stable attarctive limit cycle, whereas, the time-reversed version has a locally
% asymptotically stabel equilibrium point in the origin, which has a limit domain of
% attraction bounded by the unstable limit cycle.

% Generate symbolic state variables
syms t x1 x2 real
x = [x1;x2];

phi = [
     x1
     x2
  x1*x2
x1^2*x2
    ];

e = 1;
M = [
    0 -1 0 0 
    1 -e 0 e
    ];

f = M*phi;

%% Visualize the phase portrait of the system

scale = max([max(1,e/2),min(e,7)/1.7,min(e,1.5)/1.2]);
Visualize_2D_phase_plot2(-f,x,'XLim',[-3,3],'YLim',[-3,3]*scale,'FigNr',12)
Visualize_2D_phase_plot2(f,x,'XLim',[-3,3],'YLim',[-3,3]*scale,'FigNr',14)

%% Visualize the time evolution of the original Van der Pol oscillator system

f_ode = matlabFunction(-f,'vars',{t,x});
[tt,xx] = ode45(f_ode,[0,70],[1e-5;0]);

fig = figure(15);
Tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile; hold on, grid on, box on;
plot(tt,xx(:,1)), 
ylabel('$x_1$','Interpreter','latex','FontSize',18);
ax2 = nexttile; hold on, grid on, box on;
plot(tt,xx(:,2))
ylabel('$x_2$','Interpreter','latex','FontSize',18);
xlabel('$t$','Interpreter','latex','FontSize',18);

Link = linkprop([ax1,ax2],'XLim');

%% Derivative of phi

dot_phi = expand(jacobian(phi,x)*f);

phid = [
       x1
       x2
     x1^2
    x1*x2
     x2^2
     x1^3
  x1^2*x2
  x1*x2^2
     x2^3
     x1^4
  x1^3*x2
  x1^4*x2
    ];

Ed = find_coeff_matrix(x,phi,phid);
Ad = find_coeff_matrix(x,dot_phi,phid);

SHOULD_BE_ZERO = norm(double(simplify(  phi - Ed*phid  )))
SHOULD_BE_ZERO = norm(double(simplify(  dot_phi - Ad*phid  )))

%% Create annihilators

% N = sym(P_affine_annihilator(phi,x,'sym',1))
N = [
    x2,  0, -1,  0
     0, x1, -1,  0
     0,  0, x1, -1
     ];
[s,m] = size(N);

% Nd = sym(P_affine_annihilator(phid,x,'sym',1))
Nd = [
    x1,  0, -1,  0,  0,  0,  0,   0,   0,  0,  0,  0
    x2,  0,  0, -1,  0,  0,  0,   0,   0,  0,  0,  0
     0, x1,  0, -1,  0,  0,  0,   0,   0,  0,  0,  0
     0, x2,  0,  0, -1,  0,  0,   0,   0,  0,  0,  0
     0,  0, x1,  0,  0, -1,  0,   0,   0,  0,  0,  0
     0,  0, x2,  0,  0,  0, -1,   0,   0,  0,  0,  0
     0,  0,  0, x1,  0,  0, -1,   0,   0,  0,  0,  0
     0,  0,  0, x2,  0,  0,  0,  -1,   0,  0,  0,  0
     0,  0,  0,  0, x1,  0,  0,  -1,   0,  0,  0,  0
     0,  0,  0,  0, x2,  0,  0,   0,  -1,  0,  0,  0
     0,  0,  0,  0,  0, x1,  0,   0,   0, -1,  0,  0
     0,  0,  0,  0,  0, x2,  0,   0,   0,  0, -1,  0
     0,  0,  0,  0,  0,  0, x1,   0,   0,  0, -1,  0
     0,  0,  0,  0,  0,  0, x2, -x1,   0,  0,  0,  0
     0,  0,  0,  0,  0,  0,  0,  x2, -x1,  0,  0,  0
     0,  0,  0,  0,  0,  0,  0,   0,   0, x2,  0, -1
     0,  0,  0,  0,  0,  0,  0,   0,   0,  0, x1, -1
    ];
[sd,md] = size(Nd);

SHOULD_BE_ZERO = norm(double(simplify(  N * phi  )))
SHOULD_BE_ZERO = norm(double(simplify(  Nd * phid  )))

N_fh = matlabFunction(N,'vars',{x});
Nd_fh = matlabFunction(Nd,'vars',{x});

%%

% Corner points of the domain
X_v = [
   -0.2521    1.0387
   -1.4218    0.0202
   -0.5849   -1.2706
    1.0487   -0.8471
    1.3008    0.8370
    ];

P = sdpvar(m);
L = sdpvar(m,s,'full');
Ld = sdpvar(md,sd,'full');

CONS = [];

for i = 1:size(X_v,1)
    xi = X_v(i,:)';

    Ni = N_fh(xi);
    Ndi = Nd_fh(xi);

    CONS = [ CONS,
        P + L*Ni + Ni'*L' - eye(m) >= 0
        Ed'*P*Ad + Ad'*P*Ed + Ld*Ndi + Ndi'*Ld' + 0*eye(md) <= 0
        ];
end

sol = optimize(CONS)
P = double(P);

%% Save results


P_val = value(P);
V = phi'* P_val * phi;
V_fh = matlabFunction(V,'vars',x);

xorig_res = [51 51];
xorig_lim = [
    -1 1
    -1 1
    ]*2;

% Generate grid in X
[x1_num,x2_num] = meshgrid( ...
    linspace(-1.5,1.5,51), ...
    linspace(-1.5,1.5,51));

% Evaluate the Lyapunov function in the grid points
V_num = V_fh(x1_num,x2_num);

% Compute the maximal level set, which is still completely inside of X

lambda = linspace(0,1,101);

c1 = max(V_num(:));
c2 = 0;
for pair_of_corners = [ X_v , circshift(X_v,-1,1) ]'
    v1 = pair_of_corners(1:2);
    v2 = pair_of_corners(3:4);

    % 2D points along the line connecting v1 and v2
    p = v1*lambda + v2*(1-lambda);

    % Values of the Lyapunov function in these points
    Vp = V_fh(p(1,:),p(2,:));

    % Minimum and maximum level along the edges of X
    c1 = min([c1,Vp]);
    c2 = max([c2,Vp]);
end

% Minimum levels along the edges of the range of visualization
Mins = [min(V_num(1,:)),min(V_num(end,:)),min(V_num(:,1)),min(V_num(:,end))];

c3 = min(Mins);
c4 = max(Mins);

Min_Level = -c2;
Max_Level = 1.5*c4;

fig = figure(141);
Pc = getcontour(x1_num,x2_num,V_num,[1 1]*c1);

delete(fig.Children);
ax = axes('Parent',fig);
hold on, grid on, box on

Sf = surf(x1_num,x2_num,V_num);
Sf.CData(V_num > Max_Level) = Max_Level;
Sf.CData = log(Sf.CData + 1);
colormap turbo

X_circ = [
    X_v
    X_v(1,:)
    ];
plot3(X_circ(:,1),X_circ(:,2),X_circ*[0;0] + Min_Level,'o-','LineWidth',2)

Fl = fill3(Pc.contour(1,:),Pc.contour(2,:),Pc.contour(2,:)*0 + Min_Level,...
    min(pcz_get_plot_colors([],2)*2,1));
Fl.FaceAlpha = 0.1;
Fl.EdgeAlpha = 0;
plot3(Pc.contour(1,:),Pc.contour(2,:),Pc.contour(2,:)*0 + Min_Level,...
    'LineWidth',2,'Color',pcz_get_plot_colors([],2));

x0 = [
    -0.8
    -0.4
    ];
f_ode = matlabFunction(f,'vars',{sym('t'),x});
[tt,xx] = ode45(f_ode,[0,20],x0);

plot3(xx(:,1),xx(:,2),xx*[0;0] + Min_Level,'k')
plot3(xx(1,1),xx(1,2),Min_Level,'k.')
plot3(0,0,Min_Level,'k.','MarkerSize',15)

text(0,0,Min_Level,'origin','FontSize',12,'Interpreter','latex','VerticalAlignment','bottom','HorizontalAlignment','center')
text(xx(1,1),xx(1,2),Min_Level,'$x(0)$','FontSize',13,'Interpreter','latex','VerticalAlignment','top','HorizontalAlignment','left');
text(X_circ(4,1),X_circ(4,2),Min_Level,'$\mathcal X$','FontSize',15,'Interpreter','latex','VerticalAlignment','top','HorizontalAlignment','left')

axis([-1.5 1.5 -1.5 1.5 Min_Level-eps,Max_Level])
view([56 35])
axis vis3d

ax.TickLabelInterpreter = 'latex';
ax.FontSize = 14;

xlabel('$x_1$','Interpreter','latex','FontSize',18);
ylabel('$x_2$','Interpreter','latex','FontSize',18);
zlabel('$V(x)$','Interpreter','latex','FontSize',18);

% exportgraphics(gca,'/home/ppolcz/_/3_docs/30_ECC2021/fig/vdp_demonstrative.pdf','ContentType','vector')

