% Sol'kov system

syms t x y real

gamma = 2;
alpha = 1.2;
f = [
    1 - x*y^gamma
    alpha*y*(x*y^(gamma-1) - 1)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-5,5], ...
    'YLim',[-5,5], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('PredPray'));