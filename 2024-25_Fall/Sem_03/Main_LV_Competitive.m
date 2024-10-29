


%% Lotka-Volterra

syms t x y mu real
assumeAlso(mu >= 0.5)
assumeAlso(mu <= 2.5)
f = [
    x*(1 - x/2) - x*y
    y*(mu - y) - x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x,mu,0.5:0.1:2.5, ...
    'LaTeX',latexify(f,mu), ...
    'RunAfter',Exporter('Competitive_Mu'));

%% Lotka-Volterra

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

%%

syms t x y real
f = [
    3*x*(1 - x/5) - x*y
    y*(-4) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Mu'));

%%

syms t x y real
f = [
    3*x*(1 - x/10) - x*y
    4*y*(-1 + y/8) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Mu'));

%%

syms t x y real
f = [
    x*(3 - x/3) - x*y
    y*(-4 + y/3) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-2,15], ...
    'YLim',[-2,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Mu'));

%%

syms t x y real
f = [
    x*(-1 + 2*x) - x*y
    y*(-4 - 3*y) + 7*x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-0.2,2], ...
    'YLim',[-0.5,3], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Mu'));

%%

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
    'XLim',[-0.2,2], ...
    'YLim',[-0.5,3], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Mu'));

%%

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
    'XLim',[-0.5,5], ...
    'YLim',[-0.5,5], ...
    'LaTeX',latexify(f,p), ...
    'RunAfter',Exporter('Competitive_Mu'));

%% Simbiosis

syms t x y real
f = [
    1*x*(-1) + x*y
    2*y*(1 - y/2) + x*y
    ];
x = [x;y];

Visualize_2D_phase_plot3(f,x, ...
    'XLim',[-1,15], ...
    'YLim',[-1,15], ...
    'LaTeX',latexify(f), ...
    'RunAfter',Exporter('Competitive_Mu'));


%% Almost Van der Pol

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


%% Van der Pol

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


function ret = Exporter(Name)
    persistent DirCreated

    if isempty(DirCreated) || ~DirCreated
        mkdir("Results/" + Name)
        DirCreated = true;
    end

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