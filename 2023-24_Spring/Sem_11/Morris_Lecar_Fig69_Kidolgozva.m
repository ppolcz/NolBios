%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. November 14. (2023a)
%
% NDS project
% Morris–Lecar system

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

% parameters
V1 = -0.01; 
V2 = 0.15;
V4 = 0.04;
El = -0.5; 
Ek = -0.7; 
ECa = 1;
gl = 0.5; 
gk = 2; 
gCa = 0.9;
mu = 0.003;

syms t V w u real
x = [V;w;u];

V3 = 0.08-u;
I =  0.08-0.03*u;

m_inf = 1/2*(1+tanh((V-V1)/V2));
w_inf = 1/2*(1+tanh((V-V3)/V4));
lambda = 1/3*cosh((V-V3)/(2*V4));

% Differential equations
dV = I - gl*(V-El)-gk*w*(V-Ek)-gCa*m_inf*(V-ECa);
dw = lambda*(w_inf-w);
du = mu*(0.22+V);

f = [dV;dw;du];
f_ode = matlabFunction(f,'vars',{t,x});

x0 = [-0.275,0.0,0.0];
% x0 = [0,0.0,-0.03];
t0 = 0; tf = 2100; tspan = t0:0.01:tf; % interval

% solve the differential equation system
[t,x] = ode89(f_ode,tspan,x0); 


u_vals = sort([
    linspace(-0.043,0.0333,71)...
    linspace(-0.0337,-0.03368,101)
    ]);


V_vals1 = u_vals*0;
w_vals1 = u_vals*0;

E1 = [0;0]*u_vals;
E2 = [0;0]*u_vals*NaN;
E3 = [0;0]*u_vals*NaN;

V_vals2 = u_vals*NaN;
w_vals2 = u_vals*NaN;
V_vals3 = u_vals*NaN;
w_vals3 = u_vals*NaN;

J = jacobian([dV,dw],[V,w]);
J_fh = matlabFunction(J,'vars',[V,w,u]);

for i = 1:numel(u_vals)
    u_val = u_vals(i);
    
    dV_ = subs(dV,u,u_val);
    dw_ = subs(dw,u,u_val);

    sol1 = vpasolve([dV_,dw_],[V,w],[0.05,0.3]);
    V_vals1(i) = sol1.V;
    w_vals1(i) = sol1.w;

    E1(:,i) = eig(J_fh(double(sol1.V),double(sol1.w),u_val));

    sol2 = vpasolve([dV_,dw_],[V,w],[-0.2;0]);
    if ~isempty(sol2.V)
        V_vals2(i) = sol2.V;
        w_vals2(i) = sol2.w;
        E2(:,i) = eig(J_fh(double(sol2.V),double(sol2.w),u_val));
    end

    sol3 = vpasolve([dV_,dw_],[V,w],[-0.3;0]);
    if ~isempty(sol3.V)
        V_vals3(i) = sol3.V;
        w_vals3(i) = sol3.w;
        E3(:,i) = eig(J_fh(double(sol3.V),double(sol3.w),u_val));
    end

    fprintf('%d/%d\n',i,numel(u_vals))
end

[~,idx] = min(abs(real(E1(1,:))));

%%

fig = figure(125);
Tl = tiledlayout(4,2,"TileSpacing","compact","Padding","compact");

% plots
ax1 = nexttile;
plot(t,x(:,1),'k'); grid on;
xlim([t0,tf]); xlabel('t'); ylabel('V(t)');

ax2 = nexttile;
plot(t,x(:,1),'k'); grid on;
xlim([0,600]); xlabel('t'); ylabel('V(t)');

ax3 = nexttile;
plot(t,x(:,2),'k'); grid on;
xlim([t0,tf]); xlabel('t'); ylabel('w(t)');

ax4 = nexttile;
plot(t,x(:,3),'k'); grid on;
xlim([0,600]); xlabel('t'); ylabel('u(t)');

Link1 = linkprop([ax1,ax3],'XLim');
Link2 = linkprop([ax2,ax4],'XLim');

ax = nexttile([2,1]);
hold on; grid on; box on
plot3(x(:,3),x(:,2),x(:,1),'Color',Color_2,'LineWidth',1); 
% plot3(x(1,3),x(1,2),x(1,1),'.','Color',Color_2,'MarkerSize',15);
xlabel('u'); ylabel('w'); zlabel('V');
view([6,21]);

plot3(u_vals,w_vals1,V_vals1,'k','LineWidth',2);
plot3(u_vals,w_vals2,V_vals2,'k','LineWidth',2);
plot3(u_vals,w_vals3,V_vals3,'k','LineWidth',2);
plot3(u_vals(idx),w_vals1(idx),V_vals1(idx),'or')

ax = nexttile;
hold on; grid on; box on;
plot3(u_vals',real(E1(1,:))',imag(E1(1,:))');
plot3(u_vals',real(E2(1,:))',imag(E2(1,:))');
plot3(u_vals',real(E3(1,:))',imag(E3(1,:))');
plot(u_vals(idx),0,'or')
xlabel('u')
ylabel('real(eig)')
zlabel('imag(eig)')
view([0,90])
legend('S1 (eig1&2)','S2 (eig1)','S3 (eig1)')
xlim([-0.05,0.05])

ax = nexttile;
hold on; grid on; box on;
plot3(u_vals',real(E2(2,:))',imag(E2(2,:))','Color',Color_2);
plot3(u_vals',real(E3(2,:))',imag(E3(2,:))','Color',Color_3);
xlabel('u')
ylabel('real(eig)')
zlabel('imag(eig)')
view([0,90])
legend('S2 (eig2)','S3 (eig2)')
xlim([-0.05,0.05])

return

%%

clear x

u_vals = linspace(-0.045,0.035,51);

for u_val = u_vals
    
    dV_ = subs(dV,u,u_val);
    dw_ = subs(dw,u,u_val);
    
    dV_fh = matlabFunction(dV_,'vars',[V,w]);
    dw_fh = matlabFunction(dw_,'vars',[V,w]);

    Vw_ode = matlabFunction([dV_;dw_],'vars',{'t',[V;w]});
    
    w_range = [0,0.5];
    V_range = [-0.3,0.3];
    Vw_range = [V_range,w_range];
    
    fig = figure(123);
    Tl = tiledlayout(1,1,"Padding","compact");
    
    ax = nexttile;
    hold on, box on, grid on;
    fimplicit(dV_fh,Vw_range)
    fimplicit(dw_fh,Vw_range)
    
    sol1 = vpasolve([dV_,dw_],[V,w],[0.05,0.3]);
    S1 = plot(sol1.V,sol1.w,'o');
    
    sol2 = vpasolve([dV_,dw_],[V,w],[-0.2;0]);
    S2 = plot(sol2.V,sol2.w,'o');
    
    sol3 = vpasolve([dV_,dw_],[V,w],[-0.3;0]);
    S3 = plot(sol3.V,sol3.w,'o');
    
    [VV,ww] = meshgrid(linspace(V_range(1),V_range(2),50),linspace(w_range(1),w_range(2),50));
    dVV = dV_fh(VV,ww);
    dww = dw_fh(VV,ww);
    
    r = sqrt(dVV.^2 + dww.^2);
    dVV = dVV ./ r;
    dww = dww ./ r;
    
    Qv = quiver(VV,ww,dVV,dww,'Color',[1,1,1]*0.7);

    if ~isempty(sol1.V)
        [~,x] = ode89(Vw_ode,[0,100],double([sol1.V;sol1.w])+[1;1]*1e-3);
        plot(x(:,1),x(:,2));
    end
    if ~isempty(sol2.V)
        [~,x] = ode89(Vw_ode,[0,100],double([sol2.V;sol2.w])+[1;1]*1e-3);
        plot(x(:,1),x(:,2));
    end
    
    xlim(V_range);
    ylim(w_range);
    drawnow

    keyboard

end