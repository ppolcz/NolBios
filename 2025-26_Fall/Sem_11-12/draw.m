function [] = draw(r, varargin)


for k=1:(nargin-1)/2
    Name  = varargin{2*k-1};
    Value = varargin{2*k};
    switch lower(Name)
        case {'colorspec','legend'}
            testconfig.(lower(Name))=Value;
        otherwise
            testconfig.modelparams{i+1}  = Name;
            testconfig.modelparams{i+2}  = Value;
            i=i+2;
    end
end

color=parula(testconfig.colorspec(1));
color=color(testconfig.colorspec(2),:);

state_name = {'S', 'L', 'P', 'I', 'A', 'H'};

figure(1)
hold on; grid on;
for k=1:6
    subplot(6, 2, k); hold on; grid on;
    plot(r.t_ode, r.x_ode(:, k), 'LineWidth',2, 'Color',color);
    title(state_name{k})
end

subplot(6, 2, 7:12);
hold on; grid on;
plot(r.t_ode, r.u_ode, 'LineWidth',2, 'Color',color);
ylim([-0.5 3]);