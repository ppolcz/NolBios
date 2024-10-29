

syms t x1 x2 x3 real
x = [x1;x2;x3];

rho = 28;
sigma = 10;
beta = 8/3;

% rho = 24.0775;
rho = 30;

f = [
    sigma*(x2 - x1)
    x1*(rho - x3) - x2
    x1*x2 - beta*x3
    ];
J = jacobian(f,x);

sol = solve(f,x);

Eq = double([
    sol.x1'
    sol.x2'
    sol.x3'
    ]);

A1 = double(subs(J,x,Eq(:,1)));
[S1,D1] = eig(A1);
S1 = S1 * 10;
[~,Idx] = sort(diag(D1));
S1 = S1(:,Idx);
l1 = diag(D1(Idx,Idx));

A2 = double(subs(J,x,Eq(:,2)));
[S2,D2] = eig(A2);
l2 = diag(D2);
idx = find(abs(imag(l2)) < eps);
s2 = S2(:,idx) * 10;

A3 = double(subs(J,x,Eq(:,3)));
[S3,D3] = eig(A3);
l3 = diag(D3);
idx = find(abs(imag(l3)) < eps);
s3 = S3(:,idx) * 10;

eig(jacobian(f([1;2]),[x1;x2]))

%%

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

f_ode = matlabFunction(f,'vars',{t,x});

x0 = [0.1;3;0.11];
T = 1000;
[ttt,xxx] = ode89(f_ode,[0,T],x0);

tt = 0:0.01:T;
xx = interp1(ttt,xxx,tt);

fig = figure(123);
fig.Position(3:4) = [964 425];
Tl = tiledlayout(3,2,'Padding','compact','TileSpacing','tight');
ax = nexttile([3,1]); hold on, grid on, box on;
Pl3D = plot3(xx(:,1),xx(:,2),xx(:,3));
PlP = plot3(xx(1,1),xx(1,2),xx(1,3),'.','MarkerSize',30);
plot3(Eq(1,:),Eq(2,:),Eq(3,:),'or');
quiver3(Eq(1,1)+S1(1,1),Eq(2,1)+S1(2,1),Eq(3,1)+S1(3,1),-S1(1,1),-S1(2,1),-S1(3,1),"off",'Color',Color_5)
quiver3(Eq(1,1)-S1(1,1),Eq(2,1)-S1(2,1),Eq(3,1)-S1(3,1),+S1(1,1),+S1(2,1),+S1(3,1),"off",'Color',Color_5)
quiver3(Eq(1,1)+S1(1,2),Eq(2,1)+S1(2,2),Eq(3,1)+S1(3,2),-S1(1,2),-S1(2,2),-S1(3,2),"off",'Color',Color_5)
quiver3(Eq(1,1)-S1(1,2),Eq(2,1)-S1(2,2),Eq(3,1)-S1(3,2),+S1(1,2),+S1(2,2),+S1(3,2),"off",'Color',Color_5)
quiver3(Eq(1,1),Eq(2,1),Eq(3,1),S1(1,3),S1(2,3),S1(3,3),"off",'Color',Color_2)
quiver3(Eq(1,1),Eq(2,1),Eq(3,1),-S1(1,3),-S1(2,3),-S1(3,3),"off",'Color',Color_2)
quiver3(Eq(1,2),Eq(2,2),Eq(3,2),s2(1),s2(2),s2(3),"off",'Color',Color_5)
quiver3(Eq(1,3),Eq(2,3),Eq(3,3),s3(1),s3(2),s3(3),"off",'Color',Color_5)
quiver3(Eq(1,2),Eq(2,2),Eq(3,2),-s2(1),-s2(2),-s2(3),"off",'Color',Color_5)
quiver3(Eq(1,3),Eq(2,3),Eq(3,3),-s3(1),-s3(2),-s3(3),"off",'Color',Color_5)
plot3([0 0],[0 0],[-10,27],'k','LineWidth',2)
view([5.6 13]);
xlabel x1
ylabel x2
zlabel x3
axis equal

TlX = nexttile; hold on, box on, grid on;
plot(tt,xx(:,1));
PlX = xline(0,'r');

TlY = nexttile; hold on, box on, grid on;
plot(tt,xx(:,2));
PlY = xline(0,'r');

TlZ = nexttile; hold on, box on, grid on;
plot(tt,xx(:,3));
PlZ = xline(0,'r');

Link1 = linkprop([TlX,TlY,TlZ],'XLim');
Link2 = linkprop([PlX,PlY,PlZ],'Value');

return

%%

for i = 2:5:numel(tt)

    
    Pl3D.XData = xx(1:i,1);
    Pl3D.YData = xx(1:i,2);
    Pl3D.ZData = xx(1:i,3);

    PlP.XData = xx(i,1);
    PlP.YData = xx(i,2);
    PlP.ZData = xx(i,3);

    PlX.Value = tt(i);
    TlX.XLim = [-1,1]*10 + tt(i);

    drawnow
end


