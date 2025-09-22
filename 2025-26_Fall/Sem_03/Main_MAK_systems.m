%% Brusselator

syms a b real
syms x y real
s = [x;y];

% State equations
dx = a + x^2*y - b*x - x;
dy = -x^2*y + b*x;
f = [dx;dy];

Visualize_2D_phase_plot2(subs(f,a,1),s,b,linspace(1,3,51), ...
    'PositiveSystem',true, ...
    'XLim',[-1,5], ...
    'YLim',[-1,5], ...
    'PlotDirections',2, ...
    'LaTeX',latexify(f));

