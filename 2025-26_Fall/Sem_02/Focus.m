
M = [1 -1 ; 3 1];
A = M * [0 -1 ; 1 0] / M;
[S,D] = eig(A);

s = S(:,1) + S(:,2);


%%

[t,x] = ode45(@(t,x) A*x, [0,7], [1;0]);

fig = figure(123);
delete(fig.Children)
ax = axes(fig);
hold on;

plot(x(:,1),x(:,2))
quiver(0,0,s(1),s(2))

%%

M = [1 -1 ; 3 1];
A = M * [-2 -1 ; 1 0] / M;
[S,D] = eig(A);

s = S(:,1) + S(:,2);



fig = figure(123);
delete(fig.Children)
ax = axes(fig);
hold on;

for i = 1:10
    [t,x] = ode45(@(t,x) A*x, [0,7], randn(2,1));
    plot(x(:,1),x(:,2))
end