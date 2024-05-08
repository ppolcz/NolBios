

% Probability of infection
p = 0.01;

% Probability of recovery
q = 0.5;

% Grid size
n = 300;

% Kernel size
m = 17;

% Number of infected at day 0.
Nr_Inf = 10;

M = zeros(n,n);

Idx = ceil(rand(Nr_Inf,2)*(n-eps));
Idx = Idx(:,1) + Idx(:,2)*(n-1);
M(Idx) = 1;
M(1) = 2;

[x,y] = meshgrid(linspace(-1,1,m));
Kernel = (x.^2 + 4*y.^2) <= 1;

isS = @(M) M == 0;
isI = @(M) M == 1;

%%

fig = figure(11);
Tl = tiledlayout(1,1,"TileSpacing","tight","Padding","loose");
ax = nexttile;
hold on, grid on, box on;

colormap([
    COL.Color_1
    COL.Color_2
    COL.Color_3
    ])

surf([1,1;2,2],[1,2;1,2],[2,2;2,2]);
Sf = surf(M);
Sf.EdgeAlpha = 0;
view([0,90])
axis tight

Tl = title('k = 0');
drawnow
pause(0.1)

%%


for k = 1:1000

    % Betegek száma a szomszédban:
    nrI = conv2(isI(M),Kernel,'same');
    
    % Megfertozesi valoszinuseg
    pInf = 1 - (1-p).^nrI;
    
    % Veletlenszeru megfertozes
    newInfection = rand(n,n) < pInf & isS(M);
    
    % Veletlenszeru meggyogyulas
    newRecovery = rand(n,n) < q & isI(M);
    
    % Ha mar semmi valtozas sem tortent a populacioban, akkor hagyja abba
    if sum([newRecovery(:);newInfection(:);isI(M(:))]) == 0
        break
    end
    
    % Frissitjuk M-et
    M(newInfection) = 1;
    M(newRecovery) = 2;
    
    % Frissitjuk az abrat
    Sf.ZData = M;
    Sf.CData = M;
    Tl.String = sprintf('k = %d',k);
    drawnow

    pause(0.1)
end
