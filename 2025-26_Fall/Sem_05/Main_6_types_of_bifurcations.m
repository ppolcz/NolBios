%% Khalil, Section 2.7 Bifurcation
%
% In this script, I demostrate the six simple bifurcation types presented in Khalil's book
% in Figure 2.28: Bifurcation diagrams
% 
% The book: Khalil, Nonlinear systems (2002)
%

%% Khalil, Section 2.7, Eq. pg. 69
%  Figure 2.28(a) -- Saddle-node bifurcation

syms t x y mu real

f = [
    mu - x^2
    -y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,linspace(1,-1,9), ...
    'XLim',[-2,2], ...
    'YLim',[-2,2], ...
    ... 'LaTeX',latexify(f,mu), ...
    'RunAfter',@(varargin) keyboard),

% Visualize_2D_phase_plot2(f,x,mu,0, ...
%     'XLim',[-5,5], ...
%     'YLim',[-5,5], ...
%     'LaTeX',latexify(f,mu), ...
%     'RunAfter',@(varargin) keyboard),
%     % 'RunAfter',Exporter('Saddle-node_bifucation'));

%% Khalil, Section 2.7, Eq. on pg. 70
%  Figure 2.28(b) -- Transcritical bifurcation

syms t x y mu real

f = [
    mu*x - x^2   % This is a Lotka-Volterra model
    -y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,linspace(-1,1,9), ...
    'XLim',[-2,2], ...
    'YLim',[-2,2], ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',@(varargin) keyboard),

% Visualize_2D_phase_plot2(f,x,mu,0, ...
%     'XLim',[-5,5], ...
%     'YLim',[-5,5], ...
%     'LaTeX',latexify(f,mu), ...
%     'RunAfter',Exporter('Transcritical_bifurcation'));

%% Khalil, Section 2.7, Eq. on pg. 72
%  Figure 2.28(c) -- Supercritical pitchfork bifurcation

syms t x y mu real

f = [
    mu*x - x^3
    -y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,linspace(-1,1,9), ...
    'XLim',[-5,5], ...
    'YLim',[-5,5], ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Supercritical_pitchfork_bifurcation'));

%% Khalil, Section 2.7, Eq. on pg. 72
%  Figure 2.28(d) -- Subcritical pitchfork bifurcation

syms t x y mu real

f = [
    -mu*x + x^3 % The first term is intentionally opposed
    -y
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,linspace(-0.25,1,6), ...
    'XLim',[-2,2], ...
    'YLim',[-2,2], ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Supercritical_pitchfork_bifurcation'));

%% Khalil, Section 2.7
%  Figure 2.28(e) -- Supercritical Hopf bifurcation (Featuring Solkov's system)

syms t x y mu real

gamma = 2;
f = [
    1 - x*y^gamma
    mu*y*(x*y^(gamma-1) - 1)
    ];
x = [x;y];

Visualize_2D_phase_plot2(f,x,mu,linspace(0.9,1.4,30), ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',@(varargin) keyboard),
    ... 'RunAfter',Exporter('Selkov_system') ...
    % );


%% Khalil, Section 2.7
%  Figure 2.28(f) -- Subcritical Hopf bifurcation (time-inverted Schnakenberg's model)

syms t x y mu real

a = 1/8;
f = [
    -(a - x + x^2*y)
    -(mu - x^2*y)
    ];
xy = [x;y];

Visualize_2D_phase_plot2(f,xy,mu,linspace(0.5,0,30), ...
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',@(varargin) keyboard),
    ... 'RunAfter',Exporter('Selkov_system') ...
    % );


%% Khalil, Section 2.7, Eq. on pg. 75
%  Figure 2.30 -- Homoclinic bifurcation (global bifurcation)

syms t mu real
x1 = sym('x');
x2 = sym('y');

f = [
    x2
    mu*x2 + x1 - x1^2 + x1*x2
    ];
x = [x1;x2];

Visualize_2D_phase_plot2(f,x,mu,sort([linspace(-0.95,-0.8,30),-0.8645,-0.8645]), ...
    'XLim',[-1,2], ...
    'YLim',[-1.5,1.5], ...
    'LaTeX',latexify(f,mu) ...
    ,'RunAfter',@(varargin) keyboard ...
    ... 'RunAfter',Exporter('Selkov_system') ...
    );

%%

