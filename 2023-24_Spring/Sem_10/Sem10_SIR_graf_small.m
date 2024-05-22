%%

Ndays = 150;
simnum = 20;

%%

% Graf csucsainak szama
nV = 200;

% Minden egyes újonnan hozzáadott csúcs éleinek száma (további csúcsok hozzáadásával ez
% természetesen még nőhet).
nE_new = 3;

% Kezdeti véletlen gráf nagysága:
nV0 = 4;
nE0 = 3;
[G,A] = BAmodel(nV,nV0,nE0,nE_new);

% Barabasi-Albert graf vegso eleinek szama
nE = height(G.Edges);

if false
    [G,A] = BDmodel(nV,nE);
end

if true
    [G,A] = ERmodel(nV,nE);
end

%%

% Megjelenites (lassuk a strukturat)
[Pl1,~,ax1] = Visualize_Graph(G,1);
Pl1.EdgeAlpha = 0.5;
    
% Megjelenites (lassuk a terjedest)
[Pl2,~,ax2] = Visualize_Graph(G,2);
Pl2.EdgeAlpha = 0.1;
colormap(ax2,[
    COL.Color_1    % S
    COL.Color_Red  % I
    COL.Color_5    % R
    ])
    
fig = figure(4);
tiledlayout(1,1,"Padding","compact");
ax3 = nexttile;
hold on, grid on, box on;
title(sprintf('Fertozottek szama $(n_V = %d,n_E = %d)$',nV,nE),'Interpreter','latex','FontSize',14)

%%
% SIR modell véletlen gráfon:
% paraméterek:
p = 0.2; % 0.4, 1;
q = 0.05; % 0.2, 1;
w = 0.00001; 
perc_I_at_t0 = 0.1;
immunloss_rate = 0.2; % 1, 0.2, 0

isS = @(x) x == 0;
isI = @(x) x == 1;
isR = @(x) x == 2;

for simNr = 1:simnum
    disp(simNr)

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
        Pl2.NodeCData = x;
        drawnow
    end

    PlI = plot(ax3,I,'Color',Color,'DisplayName',Name);
    if simNr > 1
        PlI.HandleVisibility = 'off';
    end
    drawnow
end

Leg = legend('FontSize',13);