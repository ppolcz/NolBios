%  
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. October 05. (2023a)
%
function [position,isterminal,direction] = ...
    hp_ode_terminal_event_ball(t,x,P,r,args)
arguments
    t,x
    P = [0;0];
    r = 10;

    % leaving: 1
    % entering: -1
    % both directions: 0
    args.Direction {mustBeMember(args.Direction,[1,-1,0])} = 1
end

% Expression which are tested whether they approach zero.
position = norm(x - P) - r;

% Selector (toggle switch)
isterminal = 1;

% Direction of the zero crossing
direction = args.Direction;

%{

    termevent = @(t,x) hp_ode_terminal_event_ball(t,x,[1;1],0.01,'Direction',1)

%}

