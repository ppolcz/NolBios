%%

syms t x y real

mu = 1.5;

f = [
    x*(1 - x/2) - x*y
    y*(mu - y) - x*y
    ];
f_ode = matlabFunction(f,'vars',{t,[x;y]});
f1_fh = matlabFunction(f(1),'vars',[x,y]);
f2_fh = matlabFunction(f(2),'vars',[x,y]);

% Jacobi matrix szimbolikusan
J = jacobian(f,[x;y]);
% Letrehozok egy fuggvenyt a szimbolikus matrix erteku kifejezesbol
J_fh = matlabFunction(J,'vars',{[x;y]});

sol = solve(f,[x,y]);

Eq = double([ sol.x , sol.y ])';

nEq = size(Eq,2);


fig = figure(12);
delete(fig.Children)
ax = axes(fig);
hold on; grid on; box on;

plot(Eq(1,:),Eq(2,:),'.','MarkerSize',25);

dist = 0.01;

XLim = [0,3];
YLim = [0,3];

for i = 1:nEq

    eq = Eq(:,i);

    [S,lambda] = eig(J_fh(eq),'vector');

    re1 = real(lambda(1));
    re2 = real(lambda(2));
    im1 = imag(lambda(1));

    
    % Vektormezo
    [x,y] = meshgrid(linspace(XLim(1),XLim(2),101),linspace(YLim(1),YLim(2),101));
    
    f1 = f1_fh(x,y);
    f2 = f2_fh(x,y);

    r = sqrt(f1.^2 + f2.^2);
    f1 = f1 ./ r;
    f2 = f2 ./ r;

    % quiver(x,y,f1,f2)
    streamslice(x,y,f1,f2)
    
    if re1*re2 < 0
        % Nyeregpont eset

        if re1 < 0
            % Ha az elso s.e. a negativ, akkor csere
            [re1,re2] = deal(re2,re1);
            S = S(:,[2,1]);
        end

        % Stabil irany

        x0 = eq + S(:,2)*dist;
        [t,x] = ode45(f_ode,[0,-100],x0);
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);

        x0 = eq - S(:,2)*dist;
        [t,x] = ode45(f_ode,[0,-100],x0);
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);   

        % Instabil irany

        x0 = eq + S(:,1)*dist;
        [t,x] = ode45(f_ode,[0,100],x0);
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);

        x0 = eq - S(:,1)*dist;
        [t,x] = ode45(f_ode,[0,100],x0);
        plot(x(:,1),x(:,2),'r','LineWidth',1.5);
    end

end

xlim(XLim);
ylim(YLim);
