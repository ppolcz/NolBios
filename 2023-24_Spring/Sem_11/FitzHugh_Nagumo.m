
% I: parameter, x = [V;w]: allapot
syms t V w I real
x = [V;w];

c = 0.08;
b = 0.7;
a = 0.8;

dV = V - V^3/3 - w + I;
dw = c*(V + b - a*w);

% ODE-hivhato fuggveny
f = matlabFunction([dV;dw],'vars',{t,x,I});

% Racspontokban kiertekelheto fuggvenyek
dV_fh = matlabFunction(dV,'vars',{V,w,I});
dw_fh = matlabFunction(dw,'vars',{V,w,I});

% Racspontok letrehozasa
Lims = [-3,3,-3,3];
[VV,ww] = meshgrid( ...
    linspace(Lims(1),Lims(2),101), ...
    linspace(Lims(3),Lims(4),101));

T = 200;
x0 = [0.1;0];

fig = figure(123123);
fig.Position(3:4) = [1209 552];
i = 1;
for I_val = 0:0.01:1.7

    Tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile; hold on; grid on; box on;
    [tt,xx] = ode45(@(t,x) f(t,x,I_val),[0,T],x0);
    plot(tt,xx)
    
    nexttile; hold on; grid on; box on;
    fV = dV_fh(VV,ww,I_val);
    fw = dw_fh(VV,ww,I_val);

    fimplicit(@(V,w) dV_fh(V,w,I_val),Lims,'r','LineWidth',2);
    fimplicit(@(V,w) dw_fh(V,w,I_val),Lims,'r','LineWidth',2);
    streamslice(VV,ww,fV,fw);
    plot(xx(:,1),xx(:,2),'k','LineWidth',3);
    axis(Lims);

    drawnow
    % pause(0.1)
    fname = sprintf('Results/Fig%03d.jpg',i);
    exportgraphics(fig,fname)
    i = i + 1;
end