%% Khalil, Section 2.7 Bifurcation
%
% In this script, I demostrate the motion of a few particles in a specific vector field.
% This vector field has an attractive homoclinic orbit.

%% 

syms t r phi x1 x2 real

dr = r*(1-r);
dphi = sin(phi/2);

x1_expr = r*cos(phi);
x2_expr = r*sin(phi);

r_expr = sqrt(x1^2 + x2^2);
phi_expr = atan2(x2,x1) + pi;

dx = jacobian([x1_expr;x2_expr],[r;phi]) * [dr;dphi];
dx = subs(dx,[x1_expr;x2_expr],[x1;x2]);
dx = subs(dx,r,r_expr);
f = simplify(subs(dx,phi,phi_expr));

f_fh = matlabFunction(f,'vars',{t,[x1;x2]});
f1_fh = matlabFunction(f(1));          % to plot vector field
f2_fh = matlabFunction(f(2));          % to plot vector field

% Compute the values of the vector field in a given number of grid points.
[xx1,xx2] = meshgrid(linspace(-2,2,201));
f1_val = f1_fh(xx1,xx2);
f2_val = f2_fh(xx1,xx2);

% Normalize the vectors of the vector field.
r = sqrt(f1_val.^2 + f2_val.^2);
f1_val = f1_val ./ r;
f2_val = f2_val ./ r;


T = 30;
tlin = 0:0.1:T;

[tt,tr1] = ode45(f_fh,[0,T],[-1.9;-1.9]);
xlin1 = interp1(tt,tr1,tlin);

[tt,tr2] = ode45(f_fh,[0,T],[-0.1;-0.1]);
xlin2 = interp1(tt,tr2,tlin);

alpha = -pi+0.01;
[tt,tr3] = ode45(f_fh,[0,T],[cos(alpha);sin(alpha)]);
xlin3 = interp1(tt,tr3,tlin);


%%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

fig = uifigure;
g = uigridlayout(fig);
g.RowHeight = {'1x','fit'};
g.ColumnWidth = {'1x'};

ax = uiaxes(g);
streamslice(ax,xx1,xx2,f1_val,f2_val);
hold(ax,'on');
Pl3 = plot(ax,tr3(:,1),tr3(:,2),'LineWidth',2,'Color',Color_5);
Pl3 = plot(ax,tr3(1,1),tr3(1,2),'.','MarkerSize',30,'Color',Pl3.Color);
Pl1 = plot(ax,tr1(:,1),tr1(:,2),'LineWidth',2,'Color',Color_2);
Pl1 = plot(ax,tr1(1,1),tr1(1,2),'.','MarkerSize',30,'Color',Pl1.Color);
Pl2 = plot(ax,tr2(:,1),tr2(:,2),'LineWidth',2,'Color',Color_7);
Pl2 = plot(ax,tr2(1,1),tr2(1,2),'.','MarkerSize',30,'Color',Pl2.Color);

sld = uislider(g, ...
    "Limits",[0 T], ...
    "Value",0);

sld.ValueChangingFcn = @(src,event) updateRange(src,event,Pl1,Pl2,Pl3,tlin,xlin1,xlin2,xlin3);

function updateRange(src,event,Pl1,Pl2,Pl3,tlin,xlin1,xlin2,xlin3)
    x = interp1(tlin,xlin1,event.Value);
    Pl1.XData = x(1);
    Pl1.YData = x(2);
    x = interp1(tlin,xlin2,event.Value);
    Pl2.XData = x(1);
    Pl2.YData = x(2);
    x = interp1(tlin,xlin3,event.Value);
    Pl3.XData = x(1);
    Pl3.YData = x(2);
end
