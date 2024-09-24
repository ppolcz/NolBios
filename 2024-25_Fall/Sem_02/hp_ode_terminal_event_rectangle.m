%  
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2023. September 28. (2023a)
%
function [position,isterminal,direction] = ...
    hp_ode_terminal_event_rectangle(~,x,x1lim,x2lim)
arguments
    ~,x
    x1lim = [-10 10]; % Default region
    x2lim = [-10 10]; % Default region
end

% Expression which are tested whether they approach zero.
position = [
    x(1) - x1lim(1)
    x(1) - x1lim(2)
    x(2) - x2lim(1)
    x(2) - x2lim(2)
    ];

% Selector (toggle switch)
isterminal = [
    1
    1
    1
    1
    ];

% -1: zero crossing only from the negative to positive direction
% +1: zero crossing only from the positive to negative direction
%  0: zero crossing in both directions, i.e., zero can be approached from either directions
direction = [
    0
    0
    0
    0
    ];
