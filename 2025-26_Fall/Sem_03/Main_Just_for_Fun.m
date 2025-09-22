%% Competitive with Allee effect

syms t x y real
f = [
    10*x*(1 - x/8)*(x - 2) - x*y
    6*y*(1 - y/10)*(y - 2) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'PositiveSystem',true, ...
    'XLim',[-1,12], ...
    'YLim',[-1,12], ...
    'PlotDirections',1, ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_with_Allee'));

%% Homoclinic orbit
% This is only because the teacher was in a good mood when he prepared this stuff.

syms t x y real

mu = 1;
f = [
    y
    -x + mu*(1 - x^2)
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-3,3], ...
    'YLim',[-3,3], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Van_der_Pol_almost'));

%% Limit cycle, Van de Pol
% This is only because the teacher was in a good mood when he prepared this stuff.

syms t x y real

mu = 1;
f = [
    y
    -x + mu*(1 - x^2)*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x, ...
    'XLim',[-3,3], ...
    'YLim',[-3,3], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Van_der_Pol'));

%% Van der Pol (first variation)
% This is only because the teacher was in a good mood when he prepared this stuff.

syms t x y real

mu = 1;
f = [
    y
    -x + mu*(1 - x^2)*y
    ];
x = [x;y];
g = jacobian(f,x)*x;

Visualize_2D_phase_plot2(g,x, ...
    'XLim',[-3,3], ...
    'YLim',[-3,3], ...
    'LaTeX',latexify(g), ...
    'RunAfter',Exporter('Van_der_Pol_first_variation'));


%% Van der Pol (multiple variations)
% This is only because the teacher was in a good mood when he prepared this stuff.

syms t x y real

mu = 1;
f = [
    y
    -x + mu*(1 - x^2)*y
    ];
x = [x;y];

for d = 1:10
    Visualize_2D_phase_plot2(f,x, ...
        'XLim',[-3,3], ...
        'YLim',[-3,3], ...
        'LaTeX',latexify(f), ...
        'RunAfter',Exporter('Van_der_Pol_variations'));
    drawnow
    keyboard
    f = jacobian(f,x)*x;
end

%% Schnakenberg's system -- Hopf bifurcation

syms a mu x y real

assumeAlso(a > 0)
assumeAlso(mu > 0)
assumeAlso(mu < 1)

f = [
    a - x + x^2*y
    mu - x^2*y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,[a;mu],[1/8,1/8;0.1,0.2], ...
    'XLim',[0,1], ...
    'YLim',[0,3], ...
    'PlotDirections',1, ...
    'LaTeX',latexify(f,[a;mu]), ...
    'RunAfter',Exporter('Hopf_Schnakenberg'));

