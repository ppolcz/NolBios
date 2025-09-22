%% Lotka-Volterra a classical predator-prey

syms t x y real
f = [
    4*x - x*y
    -y + 0.2*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'PositiveSystem',true, ...
    'XLim',[-2,15], ...
    'YLim',[-1,8], ...
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
    'PlotDirections',2, ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_FiniteCapacity_GrowingPred'));

%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Finite and LOW carrying capacity for prey
% Predator needs prey only below a given threshold
% Coexistence with oscillating transient

syms t x y real
f = [
    x*(1 - x/5) - x*y
    -y*(1 - y/8) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-0.5,9], ...
    'YLim',[-0.5,9], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_Sc0'));


%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Finite and a CRITICAL carrying capacity for prey
% Predator needs prey only below a given threshold
% Oscillations around equilibrium

syms t x y real
f = [
    x*(1 - x/8) - x*y
    -y*(1 - y/8) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-0.5,9], ...
    'YLim',[-0.5,9], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_Sc0'));


%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Predator needs prey only below a given threshold
% Finite but HIGH carrying capacity for prey ==> opportunistic suppression
% Economic giants are the winners

syms t x y real
f = [
    3*x*(1 - x/10) - x*y
    -4*y*(1 - y/8) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-1,12], ...
    'YLim',[-1,12], ...
    'PlotDirections',2, ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_Sc1'));

%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Predator needs prey only below a given threshold
% Gradually increasing the carrying capacity of the prey
% From coexistence to opportunistic suppression

syms t x y K real
f = [
    x*(1 - x/K) - x*y
    -y*(1 - y/3) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,K,linspace(2,4,51), ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,K), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_Bifurcation'));

%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Predator needs prey only below a given threshold
% Gradually decreasing the predation coefficient
% Opportunistic suppression of poor to exclusion of the rich

syms t x y a real
f = [
    x*(1 - x/4) - a*x*y
    -y*(1 - y/3) + a*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,a,linspace(1,0.1,51), ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,a), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_Bifurcation2'));

%% Lotka-Volterra predator-prey
% Predator: economic giants
% Prey: common citizens
% Predator needs prey only below a given threshold
% Very low predation coefficient
% Commom citizens can be the winners

syms t x y a real
f = [
    x*(1 - x/4) - a*x*y
    -y*(1 - y/3) + a*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,a,0.1, ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,a), ...
    'RunAfter',Exporter('Economic_Giants_vs_Common_citizens_SmallDep'));

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
    'PositiveSystem',true, ...
    'XLim',[-4,20], ...
    'YLim',[-3,15], ...
    'PlotDirections',2, ...
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
    'PositiveSystem',true, ...
    'XLim',[-1,60], ...
    'YLim',[-1,60], ...
    'PlotDirections',10, ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray_DyingPray_LimitedStomach'));

