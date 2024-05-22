%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. November 14. (2023a)
%
% NDS project
% FitzHugh-Nagumo system

Color_1 = [0 0.4470 0.7410];
Color_2 = [0.8500 0.3250 0.0980];
Color_3 = [0.9290 0.6940 0.1250];
Color_4 = [0.4940 0.1840 0.5560];
Color_5 = [0.4660 0.6740 0.1880];
Color_6 = [0.3010 0.7450 0.9330];
Color_7 = [0.6350 0.0780 0.1840];

syms t V w u real

a = -u;
b = 1.3;
c = -0.32;
d = 0.05;
e = 1;
mu = -0.0001;

dV = V - V^3/3 - w + u;
dw = 0.08*(V + 0.7 - 0.8*w);
du = mu*(0.7 + V);

f = [dV;dw;du];
f_ode = matlabFunction(f,'vars',{t,[V;w;u]});

Tp = 141.76 * 10;
To = 42.4;
x0 = [-0.3122,-0.3530,-0.7640];
x0 = [-1.59056,-0.2226,-0.982981];
x0 = [-1.08846,-0.669485,0.35];

% solve the differential equation system
[t,x] = ode89(f_ode,[0,3*Tp+To],x0); 

u_val = -1;
% Find initial guess for the stationary points
%{
    dV_fh = matlabFunction(subs(dV,u,u_val),'vars',[V,w]);
    dw_fh = matlabFunction(subs(dw,u,u_val),'vars',[V,w]);

    w_range = [-2,2];
    V_range = [-2,2];
    Vw_range = [V_range,w_range];
    
    fig = figure(123);
    Tl = tiledlayout(1,1,"Padding","compact");
    
    ax = nexttile;
    hold on, box on, grid on;
    fimplicit(dV_fh,Vw_range)
    fimplicit(dw_fh,Vw_range)
%}

S = {
    [0.55;0.48]
    [-0.4;-0.38]
    [-0.75;-0.6]
    };

getsol = @(sol) double([sol.V;sol.w]);
mysolve = @(dV,dw,x0) vpasolve([dV;dw],[V;w],x0);

u_vals = sort([
    linspace(-1,-0.75,71) ...
    % linspace(-0.0337,-0.03368,101)
    ]);

SS = repmat({[0;0]*u_vals*NaN},[3,1]);
EE = repmat({[0;0]*u_vals*NaN},[3,1]);

J = jacobian([dV,dw],[V,w]);
J_fh = matlabFunction(J,'vars',{[V;w],u});

% for i = 1:numel(u_vals)
%     u_val = u_vals(i);
% 
%     dV_ = subs(dV,u,u_val);
%     dw_ = subs(dw,u,u_val);
% 
%     for s = 1:3
%         sol = mysolve(dV_,dw_,S{s});
% 
%         if ~isempty(sol.V)
%             S{s} = getsol(sol);
%             SS{s}(:,i) = S{s};
%             EE{s}(:,i) = eig(J_fh(S{s},u_val));
%         end
%     end
% 
%     fprintf('%d/%d\n',i,numel(u_vals));
% end

%%

fig = figure(125);
Tl = tiledlayout(4,2,"TileSpacing","compact","Padding","compact","TileIndexing","columnmajor");

% plots
axV = nexttile;
plot(t,x(:,1),'k'); grid on;
xlim([t(1),t(end)]); xlabel('t'); ylabel('V(t)');

axw = nexttile;
plot(t,x(:,2),'k'); grid on;
xlim([t(1),t(end)]); xlabel('t'); ylabel('w(t)');

axu = nexttile;
plot(t,x(:,3),'k'); grid on;
xlim([t(1),t(end)]); xlabel('t'); ylabel('u(t)');

Link1 = linkprop([axV,axw,axu],'XLim');

ax = nexttile([2,1]);
hold on; grid on; box on
plot3(x(:,3),x(:,2),x(:,1),'Color',Color_2,'LineWidth',1); 
% plot3(x(1,3),x(1,2),x(1,1),'.','Color',Color_2,'MarkerSize',15);
xlabel('u'); ylabel('w'); zlabel('V');
view([6,21]);

for s = 1:3
    plot3(u_vals,SS{s}(2,:),SS{s}(1,:),'k','LineWidth',2);
end

ax = nexttile;
hold on; grid on; box on;
for s = 1:3
    plot3(u_vals,real(EE{s}(1,:)),imag(EE{s}(1,:)));
end
xlabel('u')
ylabel('real(eig)')
zlabel('imag(eig)')
legend('S1 (eig1&2)','S2 (eig1)','S3 (eig1&2)')

ax = nexttile;
hold on; grid on; box on;
plot3(u_vals,real(EE{2}(2,:)),imag(EE{2}(2,:)),'Color',Color_2);
xlabel('u')
ylabel('real(eig)')
zlabel('imag(eig)')
legend('S2 (eig2)')


return

%%

clear x

u_vals = linspace(-1,-0.75,51);

w_range = [-4,4];
V_range = [-4,4];

for u_val = u_vals
    
    dV_ = subs(dV,u,u_val);
    dw_ = subs(dw,u,u_val);
    
    dV_fh = matlabFunction(dV_,'vars',[V,w]);
    dw_fh = matlabFunction(dw_,'vars',[V,w]);

    Vw_ode = matlabFunction([dV_;dw_],'vars',{'t',[V;w]});
    
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