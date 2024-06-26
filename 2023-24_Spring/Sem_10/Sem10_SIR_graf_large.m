%%

% Graf csucsainak szama
nV = 10000;

% Minden egyes újonnan hozzáadott csúcs éleinek száma (további csúcsok hozzáadásával ez
% természetesen még nőhet).
nE_new = 2;

% Kezdeti véletlen gráf nagysága:
nV0 = 4;
nE0 = 3;
[G1,A1] = BAmodel(nV,nV0,nE0,nE_new);

% Barabasi-Albert graf vegso eleinek szama
nE = height(G1.Edges);

%% Egy teljesen rendezetlen véletlen gráf

A2 = sparse(nV,nV);
idx = find(tril(A2+1,-1));
sigma = randperm(numel(idx),nE);
A2(idx(sigma)) = 1;
A2 = A2 + A2' - diag(A2);
G2 = graph(A2);

%% Blokk diagonális véletlen gráf

A3 = sparse(nV,nV);
idx = find(tril(A3+1,-1));
sigma = randperm(numel(idx),nE);
A3(idx(sigma)) = 1;
A3 = A3 + A3' - diag(A3);
G3 = graph(A3);

%%

fig = figure(4);
tiledlayout(1,1,"Padding","compact");
ax2 = nexttile;
hold on, grid on, box on;
title(sprintf('Fertozottek szama $(n_V = %d,n_E = %d)$',nV,nE),'Interpreter','latex','FontSize',14)

Ndays = 100;
simnum = 20;

% SIR modell véletlen gráfon:
% paraméterek:
p = 0.1; % 0.1, 0.4, 1;
q = 0.1; % 0.2, 1;
w = 0.0027; 
perc_I_at_t0 = 0.005;
immunloss_rate = 0.2; % 1, 0.2, 0

syms S I R beta t real
dS = -beta*S*I/nV + w*R;
dI = beta*S*I/nV - q*I;
dR = q*I - w*R;
f = matlabFunction([dS;dI;dR],'vars',{t,[S;I;R],beta});

x0 = [
    1-perc_I_at_t0
    perc_I_at_t0
    0
    ]*nV;

beta = 1-(1-p).^mean(degree(G1));
[t,x] = ode45(@(t,x) f(t,x,beta),[0,Ndays],x0);
plot(t,x(:,2),'Color',COL.Color_Black,'LineWidth',2,'DisplayName','ODE solution')


for Selected_Graph = [1,2,3] 
    
    switch Selected_Graph
        case 1
            G = G1;
            A = A1;
            Color = COL.Color_1;
            Name = 'Barabasi-Albert modell';
        case 2
            G = G2;
            A = A2;
            Color = COL.Color_2;
            Name = 'Erdos-Renyi modell';
        case 3
            G = G3;
            A = A3;
            Color = COL.Color_3;
            Name = 'Blokk-diagonalis modell';
    end
    
    %%
    
    if nV <= 1000
        % Megjelenites
        [Pl,~,ax] = Visualize_Graph(G);
        Pl.EdgeAlpha = 0.5;
    else
        Pl = [];
    end
        
    %%
    
    isS = @(x) x == 0;
    isI = @(x) x == 1;
    isR = @(x) x == 2;
    
    % colormap(ax,[
    %     COL.Color_1    % S
    %     COL.Color_Red  % I
    %     COL.Color_5    % R
    %     ])
    % Pl.EdgeAlpha = 0.1;
    
    for simNr = 1:simnum
        % disp(simNr)
    
        % itt tárolom a csúcsok állapotát: 
        % 0 egészséges, de fogékony, 
        % 1 fertõzött, 
        % 2 felgyógyult immunis
        x = zeros(nV,1);
        S = zeros(Ndays+1,1);
        I = zeros(Ndays+1,1);
        
        x(rand(size(x)) <= perc_I_at_t0) = 1;
        I(1) = sum(isI(x));
        S(1) = sum(isS(x));
        for i = 1:Ndays % a járvány lépései
        
            % Betegek száma a szomszédban:
            nrI = A * isI(x);
            
            % Megfertozesi valoszinuseg
            pInf = 1 - (1-p).^nrI;
            
            % Veletlenszeru megfertozes
            newInfection = rand(size(x)) < pInf & isS(x);
            
            % Veletlenszeru meggyogyulas
            newRecovery = rand(size(x)) < q & isI(x);
            
            % Veletlenszeru meggyogyulas
            newImmunityLoss = rand(size(x)) < w & isR(x);
            
            % Ha mar semmi valtozas sem tortent a populacioban, akkor hagyja abba
            if sum(isI(x)) == 0
                break
            end
            
            % Frissitjuk a csucsok erteket
            x(newInfection) = 1;
            x(newRecovery) = 2;
            x(newImmunityLoss) = 0;
    
            I(i+1) = sum(isI(x));
            S(i+1) = sum(isS(x));
            
            % Frissitem a graf szinertekeit
            % Pl.NodeCData = x;
            % drawnow
        end
    
        PlI = plot(I,'Color',Color,'DisplayName',Name);
        if simNr > 1
            PlI.HandleVisibility = 'off';
        end
        drawnow
    end
end

Leg = legend('FontSize',13);