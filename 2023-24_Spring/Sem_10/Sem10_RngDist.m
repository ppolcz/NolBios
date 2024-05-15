%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. May 15. (2023a)
%


rng(1)

%% Random sample generator for the random input variable

N = 10;

k = 1:N;
P = k.^(-3);
P = P / sum(P);

F = [0 cumsum(P)];
p = linspace(0,1,N+1);

MM = 100000;
XX = rand(1,MM);
YY = interp1(F,p,XX) * N;

M = 3;
X = rand(1,M);
Y = ceil(interp1(F,p,X) * N);

fig = figure(1);
delete(fig.Children)
ax = axes(fig);
hold on, box on, grid on

histogram(YY,N,'Normalization','pdf')
plot(Y-0.5,randn(1,M)*0.1+0.5,'*r')
