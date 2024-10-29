

syms t x1 x2 x3 real
x = [x1;x2;x3];

b = 0.18;
f = [
    -b*x1 + sin(x2)
    -b*x2 + sin(x3)
    -b*x3 + sin(x1)
    ];

f_ode = matlabFunction(f,'vars',{t,x});

x0 = [0.1;-3;0.11];
T = 1000;
[tt,xx] = ode89(f_ode,[0,T],x0);

fig = figure(123);
delete(fig.Children);
ax = axes(fig); hold on
plot3(xx(:,1),xx(:,2),xx(:,3))