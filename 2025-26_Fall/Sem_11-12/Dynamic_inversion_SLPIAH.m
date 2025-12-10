% clear; close all; clc;

% set(0,'DefaultFigureWindowStyle','docked')

%% Model definition

% Declare symbolic variables:
% beta: transmission rate of the pathogen (the unknown external signal)
syms beta nu omega real
syms S L P I A H D R real
syms tau_L tau_P tau_I tau_A tau_H p_I p_H p_D q_A real
syms N_p real

% The state dynamics
dS = -beta*(P+I+q_A*A)*S - nu*S + omega*R;
dL = beta*(P+I+q_A*A)*S - tau_L*L;
dP = tau_L*L - tau_P*P;
dI = p_I*tau_P*P - tau_I*I;
dA = (1-p_I)*tau_P*P - tau_A*A;
dH = tau_I*p_H*I - tau_H*H;
dD = p_D*tau_H*H;
dR = tau_I*(1-p_H)*I + tau_A*A + (1-p_D)*tau_H*H + nu*S - omega*R;

u = beta;
x = [S;L;P;I;A;H;D;R];
F_sym = [dS;dL;dP;dI;dA;dH;dD;dR];

f_sym = subs(F_sym,beta,0);
g_sym = subs(F_sym - f_sym,u,1);
h = H;

% Parameter values
Par.nu = 0;        % vaccination rate
Par.omega = 1/120; % expected time of losing immunity is about 4 month (120 days)
Par.tau_L = 1/2.5; % expected length of the latent phase is 2.5 days
Par.tau_P = 1/3;   % expected length of the presymptomatic phase: 3 days
Par.tau_I = 1/4;   % expected length of the main phase of the disease (with symptoms): 4 days
Par.tau_A = 1/4;   % expected length of the main phase of the disease (WITHOUT symptoms): 4 days
Par.tau_H = 1/10;  % expected length of hospital treatment: 10 days
Par.p_I = 0.6;     % probability of showing symptoms in the main phase: 60%
Par.p_H = 0.076;   % probability that hospital treatement needed: 7.6%
Par.p_D = 0.02;    % probability of fatal outcome if hospitalized
Par.q_A = 0.75;    % relative infectiousness of asymptomatic people
Par.N_p = 9800000; % population

P = struct2table(Par);
params = sym(P.Properties.VariableNames);
values = P.Variables;

f = matlabFunction(subs(f_sym,params,values), 'Vars', {x});
g = matlabFunction(subs(g_sym,params,values), 'Vars', {x});
F = @(x,u) f(x) + g(x)*u;

%% Create control model and calculate U depending on Reference

syms u r real
y = h;

Lg = @(h) jacobian(h,x)*g_sym;
Lf = @(h) jacobian(h,x)*f_sym;

Lfh = Lf(h)
Lgh = Lg(h) % is zero
dy = Lfh + Lgh*u

Lf2h = Lf(Lfh)
LgLfh = Lg(Lfh) % is zero
ddy = Lf2h + LgLfh*u

Lf3h = Lf(Lf2h)
LgLf2h = Lg(Lf2h) % is zero
dddy = Lf3h + LgLf2h*u

Lf4h = Lf(Lf3h)
LgLf3h = Lg(Lf3h) % is NOT zero
d4y = Lf4h + LgLf3h*u

% The linearizing feedback is:
% u = -(Lf4h - v)/LgLf3h
% but first, we need to design a controller, which manipulates the
% simulator through the new signal v.

% The system from input v to output y acts like a linear system of the form
A = [0 1 0 0;
     0 0 1 0;
     0 0 0 1;
     0 0 0 0];
B = [0; 0; 0; 1];
C = [1, 0, 0, 0];

% The state of this linear system (xi) can be given as a partial
% transformation of the state (x) of the original system:
xi = [
    h
    Lfh
    Lf2h
    Lf3h
    ];

% Design a pole-placement controller for the linear system
K = place(A, B, [-3 -4 -5 -6]); % feedback gain
G = -inv(C/(A-B*K)*B);          % feedforward gain

% Control signal:
v = G*r-K*xi;

% The linearizing feedback with the control signal is:
u = -(Lf4h - v)/LgLf3h

u = matlabFunction(subs(u,params,values),'Vars',{x,r});

%% Load and fit data

% Load hospitalization data and normalize
L = load("H_All.mat");
H_All = L.H_All;
H_All = H_All/Par.N_p;

% The data is noise, therefore, we first fit spline polynomials to the data
H_pp = spline(0:length(H_All)-1,H_All);
H_measured = @(t) ppval(H_pp,t);

START = 153;
STOP = 455;

H0 = H_measured(START);

% Simulate
x0 = [0; 3*H0; 8*H0; 5*H0; 4*H0; H0;0;0];
x0(1) = 1-sum(x0(2:end));

f_ode = @(t,x) F(x,u(x,H_measured(t)));
tspan = [START,STOP];
[tt,xx] = ode89(f_ode,tspan,x0);

uu = u(xx',H_measured(tt'))';

%%

fig1 = figure(1);
delete(fig1.Children)
hold on, grid on, box on;
plot(tt-START, xx(:, 2:end), 'LineWidth',3);
xlim([0 300]);
ylim([0 10*1e-3])
xlabel('Day');
ylabel('Number of people (relative to population)')
Leg1 = legend('L','P', 'I', 'A', 'H','D','R', 'Location', 'northwest');
ax1 = gca;
title('Reconstructed state of the epidemic', ...
    'Interpreter','latex','FontSize',14)

fig2 = figure(2);
delete(fig2.Children)
hold on, grid on, box on;
plot(tt-START, xx(:, 6), 'LineWidth',3);
plot(tt-START, H_measured(tt),'LineWidth',2)

Leg2 = legend('Estimated H(t)', 'H_{official}', 'Location', 'northwest');
xlim([0 300]);
xlabel('Day');
ylabel('Number of people (relative to population)')
ax2 = gca;

fig3 = figure(3);
delete(fig3.Children)
hold on, grid on, box on;
plot(tt-START, uu, 'LineWidth',3);
title('Transmission rate of the pathogen, $\beta$: the unknown input', ...
    'Interpreter','latex','FontSize',14)
xlim([0 300]);
ylim([0 1]);
xlabel('Day');
ax3 = gca;

for ax = [ax1 , ax2 , ax3]
    ax.FontSize = 20;
    ax.TickLabelInterpreter = 'latex';
end

for Leg = [Leg1 , Leg2]
    Leg.Interpreter = 'latex';
    Leg.FontSize = 20;
end
