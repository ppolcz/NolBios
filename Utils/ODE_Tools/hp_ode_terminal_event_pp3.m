%  
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. October 05. (2023a)
%
function [position,isterminal,direction] = ...
    hp_ode_terminal_event_pp3(t,x,x1lim,x2lim,x3lim,P,r)
arguments
    t
    x (3,1)
    x1lim,x2lim,x3lim
    P (3,:)
    r = 0.1;
end

% Expression which are tested whether they approach zero.
position = [
    x(1) - x1lim(1)
    x(1) - x1lim(2)
    x(2) - x2lim(1)
    x(2) - x2lim(2)
    x(3) - x3lim(1)
    x(3) - x3lim(2)
    vecnorm(x-P,1)' - r
    ];

% Selector (toggle switch)
isterminal = [
    1
    1
    1
    1
    1
    1
    ones(size(P,2),1)
    ];

% Direction of the zero crossing
direction = [
    -1
    1
    -1
    1
    -1
    1
    -ones(size(P,2),1)
    ];
