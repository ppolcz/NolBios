%% Lotka-Volterra, competitive asymmetric

syms t x y real
f = [
    3*x*(1 - x/10) - 0.15*x*y
    y*(1 - y/10)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'PlotDirections',2, ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Asymmetric'));


%% Lotka-Volterra, competing, mutually exclusive species

syms t x y mu real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
f = [
    x*(1 - x/2) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x,mu,0.5:0.1:2.5, ...
    'PositiveSystem',true, ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Competitive_Mu'));

%% Lotka-Volterra, competitive coexistence

syms t x y mu real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
f = [
    x*(2 - 2*x) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x,mu,0.5:0.1:2.5, ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Competitive_Mu'));

%% Lotka-Volterra, competitive species
%  From exclusion to coexistence (parametrization 1.)

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

Visualize_2D_phase_plot3(f,x,p,p_vals, ...
    'PositiveSystem',true, ...
    'XLim',[-0.2,2], ...
    'YLim',[-0.5,3], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Mu'));

%% Lotka-Volterra, competitive species
%  From exclusion to coexistence (parametrization 2.)
%  Fixed equilibrium.

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

Visualize_2D_phase_plot3(f,x,p,p_vals, ...
    'PositiveSystem',true, ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Mu'));
