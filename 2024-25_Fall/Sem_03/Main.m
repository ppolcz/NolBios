%% Lotka-Volterra predator-prey

syms t x y real
f = [
    4*x - x*y
    -y + 0.2*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray'));

%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey

syms t x y real
f = [
    4*x*(1 - x/10) - x*y
    -y + 0.2*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_FiniteCapacity'));


%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey
% Predator needs prey only below a given threshold

syms t x y real
f = [
    3*x*(1 - x/10) - x*y
    -y*(1 - y/8) + 0.3*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,12], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_FiniteCapacity_GrowingPred'));

%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey
% Predator needs prey only below a given threshold

syms t x y real
f = [
    3*x*(1 - x/10) - x*y
    4*y*(-1 + y/8) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_FiniteCapacity_GrowingPred2'));

%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey
% Predator needs prey only below a given threshold

syms t x y real
f = [
    x*(3 - x/3) - x*y
    y*(-4 + y/3) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_FiniteCapacity_GrowingPred3'));


%% Lotka-Volterra predator-prey
% Preys cannot live below a given population

syms t x y real
f = [
    x*(-1 + 2*x) - x*y
    y*(-4 - 3*y) + 7*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-0.2,2], ...
    'YLim',[-0.5,3], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_DyingPray'));


%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey
% Predator cannot eat infinitely many prey

a = 10;
D = 4;
syms t x y real
f = [
    -3*x*(1 - x/3) - a*x*y/(D + x)
    -4*y*(1 + 3*y/4) + 3*a*x*y/(D + x)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-4,20], ...
    'YLim',[-3,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_DyingPray_LimitedStomach'));

%% Lotka-Volterra predator-prey
% Finite carrying capacity for prey
% Predator cannot eat infinitely many prey
% Limit cycle appears

r = 1.5;
a = 5;
b = 18;
D = 2;
c = 4;

K_Hopf = b*(a*D+c)/(a*D-c);
K = K_Hopf * (0.8 + 0.4);
K = 50;

syms t x y real
f = [
    r*x*(1 - x/K) - y*(a*x)/(b+x)
    -c*y + D*y*(a*x)/(b+x)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-6,100], ...
    'YLim',[-6,100], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_DyingPray_LimitedStomach'));

%% Lotka-Volterra competitive
% Asymmetric

syms t x y real
f = [
    3*x*(1 - x/10) - 0.15*x*y
    y*(1 - y/10)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Asymmetric'));


%% Lotka-Volterra competitive
% Symmetric, parametric

syms t x y mu real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
f = [
    x*(1 - x/2) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,0.5:0.1:2.5, ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Competitive_Exclusion_Mu'));

%% Lotka-Volterra competitive
% Symmetric, parametric

syms t x y mu real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
f = [
    2*x*(1 - x) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,0.5:0.1:2.5, ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Competitive_Coexistence_Mu'));

syms t x y mu eta real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
assumeAlso(eta >= 1)
assumeAlso(eta <= 2)
f = [
    eta*x*(1 - eta*x/2) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];
p = [mu;eta];

eta_vals = 1:0.05:2;
p_vals = [
    eta_vals*0 + 1.3
    eta_vals
    ];

Visualize_2D_phase_plot2(f,x,p,p_vals, ...
    'XLim',[-0.2,2], ...
    'YLim',[-0.5,3], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Rotating'));

%% Lotka-Volterra competitive
% Symmetric, parametric

xe = 1;
ye = 2;

syms t x y mu K real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
assumeAlso(K >= 1)
assumeAlso(K <= 5)
f = [
    x*( ye*(1-x/K) - y*(1-xe/K) )
    y*( xe*(1-y/K) - x*(1-ye/K) )
    ];
x = [x;y];
p = [mu;K];

K_vals = 1.5:0.1:4.5;
p_vals = [
    K_vals*0 + 1.3
    K_vals
    ];

Visualize_2D_phase_plot2(f,x,p,p_vals, ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Rotating0'));

%% Simbiosis

syms t x y real
f = [
    1*x*(-1) + x*y
    2*y*(1 - y/2) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-1,15], ...
    'YLim',[-1,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Simbiosis'));

%% Simbiosis

syms t x y real
f = [
    x*(5-x) + 0.7*x*y
    y*(5-y) + 0.7*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-1,20], ...
    'YLim',[-1,20], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Simbiosis2'));

%% Competitive with Allee effect

syms t x y real
f = [
    10*x*(1 - x/8)*(x - 2) - x*y
    6*y*(1 - y/10)*(y - 2) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_with_Allee'));

%% Functions

function ret = Exporter(Name)
    mkdir("Results/" + Name)
    ret = @(Fig,ind,p_val) exportgraphics(Fig,sprintf('Results/%s/%03d.jpg',Name,ind));
end

function ret = latexify(f,p)
    arguments
        f, p = []
    end

    Eq = "$\left\{\begin{array}{l} \dot x = " + latex(f(1)) + "\\ \dot y = " + latex(f(2)) + " \end{array}\right.$";

    if isempty(p)
        ret = @(~) Eq;
        return
    end

    pnames = cellfun(@(s) {[', $' latex(s) ' = ']}, num2cell(p));

    function ret = lambdaFun(p_val)
        ret = Eq + strjoin(cellfun(@(s,v) string(s) + num2str(v) + "$",pnames,num2cell(p_val)),'');
    end

    ret = @lambdaFun;
end