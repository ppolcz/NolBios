%%
% Author: Peter Polcz (ppolcz@gmail.com)  
% Created on 2023. September 28. (2023a)

A = [-1 2 ; 2 -1];
f = @(t,x) A*x + [
    0.7*x(1)^2
    -1.8*x(1)*x(2)
    ];

T = 10;

x1lim = [-3 3];
x2lim = [-3 3];
term_event = @(t,x) hp_ode_termevent_rectangle(t,x,x1lim,x2lim);
odeopts = odeset('Events',term_event);

%{
Although the following solution is formally correct, anonymous functions implemented by
the built-in function `deal` are not supported by the ODE solvers. Therefore, I created my
own dealer function in a separate file `hp_dealer_function.m`, which contains:

  function [varargout] = hp_dealer_function(varargin) 
      varargout = varargin(1:nargout); 
  end

%}

% This should work
%{
term_event = @(t,x) my_dealer_function([
    x(1) - x1lim(1)
    x(1) - x1lim(2)
    x(2) - x2lim(1)
    x(2) - x2lim(2)
    ],[ 1 1 1 1 ]',[ 0 0 0 0 ]');
odeopts = odeset('Events',term_event);
%}



x0_From = [1;-3];
x0_To = [2;-2];
lambda = linspace(0,1,500);

fig = figure(1);
delete(fig.Children)
ax = axes(fig);
hold on, grid on, box on

for x0 = x0_From*lambda + x0_To*(1-lambda)
    [t_sol,x_sol] = ode45(f,[0 T],x0,odeopts);

    if x_sol(end,2) > 0
        Color = [1 0 0];
    else
        Color = [0 1 0];
    end

    Pl = plot(x_sol(:,1),x_sol(:,2),'Color',Color);
    drawnow
end

x0_From = [-1;1];
x0_To = [-1.5;0.5];

for x0 = x0_From*lambda + x0_To*(1-lambda) ,
    [t_sol,x_sol] = ode45(f,[0 T],x0,odeopts);

    if x_sol(end,2) < 0
        Color = [1 0 0];
    else
        Color = [0 1 0];
    end

    Pl = plot(x_sol(:,1),x_sol(:,2),'Color',Color);
    drawnow
end

%%

[S,D] = eig(A);
S = S*D;
s1 = S(:,1);
s2 = S(:,2);

quiver(s1(1),s1(2),-s1(1),-s1(2),1,"LineWidth",3,'Color',[0 0 0])
quiver(0,0,s2(1),s2(2),1,"LineWidth",3,'Color',[0 0 0])

axis equal

%%

function [varargout] = my_dealer_function(varargin)
    varargout = varargin(1:nargout);
end
