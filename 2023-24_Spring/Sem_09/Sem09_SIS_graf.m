%P�lda gr�fok:
% [G] = RegularNetwork(1,0, 50, 2, 0.0); %0; 0.2; 0.7; 1
[G] = ERmodel(0,0, 50, 0, 120);
% [G] = BAmodel(1,1, 50, 3, 2, 3);

itnum = 150;
simnum = 20;

%%
% Gr�frajzol� subrutin. Elhelyezi a gr�f nv cs�cs�t egy k�rvonal ment�n
% Created by Pablo Blinder. blinderp@bgu.ac.il
% Updated by David Tekan. tekda@digitus.itk.ppke.hu
center=[0,0];
theta=linspace(0,2*pi,G.nv+1);
rho=ones(1,G.nv+1);
[X,Y] = pol2cart(theta',rho');
X=X+center(1);
Y=Y+center(2);
x=(X(1:end-1)*10)';
y=(Y(1:end-1)*10)';
%
h = nan(G.nv,1);

fig = figure(13);
fig.Position(3:4) = [1157,520];
Tl = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
ax = nexttile;
[XX,YY]=gplot(G.Adj,[x' y'],'k-');
i=~isnan(XX);
XX=XX(i);YY=YY(i);
XX=reshape(XX,2,length(XX)/2);
YY=reshape(YY,2,length(YY)/2);
hLines=line(XX,YY);
set(hLines,'color','k');

hold on;
kv=full(diag(G.Adj*G.Adj));
kvGroups=unique(setdiff(kv,0));
map=jet(max(kvGroups));

for i=1:G.nv
    h(i,:)=plot(x(i),y(i),'ko');
    text(x(i)+0.1*x(i),y(i)+0.1*y(i),num2str(i));
end

set(h,'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerSize',10);
set(gca,'Visible','Off','YDir','reverse');
colormap(map);
hc=colorbar;

set(hc,'FontSize',8,'FontW','Demi')
set(hc,'Visible','off')
set(gcf,'Color','w')

%fert�z�ttek ar�ny�nak �br�zol�sa
ax = nexttile;
hold on, grid on, box on;

immunvesztes_idopontok = round( [35,70]/100*itnum );

fert_num = nan(simnum,itnum+1); 
Pl = plot(0:itnum,fert_num);
Xl = xline(immunvesztes_idopontok,'-',repmat({'immunvesztes'},size(immunvesztes_idopontok)));
ylim([0 1]);
xlim([0 itnum]);
xlabel('l�p�s');
ylabel('fert�z�tt');

% atlag_fert = mean(fert_num(itnum/2:end)) %�tlag fert�z�tts�gi szint a j�rv�ny 2. fel�ben (tranziensek ne sz�m�tsanak)

%%
% SIR modell v�letlen gr�fon:
% param�terek:
p = 0.2; % 0.4, 1;
q = 0.05; % 0.2, 1;
w_new_variant = 0.7;
w_standard = 0.02; 
kezdetifert = 0.01;
immunloss_rate = 0.2; % 1, 0.2, 0

isS = @(csucs) csucs == 0;
isI = @(csucs) csucs == 1;
isR = @(csucs) csucs == 2;

% waitforbuttonpress;

for simNr = 1:simnum
    disp(simNr)

    % itt t�rolom a cs�csok �llapot�t: 
    % 0 eg�szs�ges, de fog�kony, 
    % 1 fert�z�tt, 
    % 2 felgy�gyult immunis
    csucs = zeros(G.nv,1); 
    
    kezd = rand(size(csucs));
    fert = find(kezd <= kezdetifert); % alapb�l ennyi fert�z�tt cs�ccsal kezd�nk
    csucs(fert) = 1;
    hfert = h(fert); % fert�z�tt cs�csok koordin�t�i
    
    egeszs = find(csucs==0); % nem fert�z�ttek
    hegeszs = h(egeszs); % eg�szs�ges cs�csok koordin�t�i
    
    fert_num(simNr,1)=length(fert)/G.nv; %fert�z�ttek ar�nya az �sszes cs�cshoz k�pest
    %a cs�csok sz�nez�se az �llapotuk szerint: fert: piros, eg�szs: z�ld, immunis: k�k
    set(hegeszs,'MarkerFaceColor',[0 1 0]); 
    set(hfert,'MarkerFaceColor',[1 0 0]);
    
    for i = 1:itnum % a j�rv�ny l�p�sei
    
        if ismember(i,immunvesztes_idopontok)
            w = w_new_variant;
        else
            w = w_standard;
        end
    
        % Betegek sz�ma a szomsz�dban:
        nrI = G.Adj * isI(csucs);
        
        % Megfertozesi valoszinuseg
        pInf = 1 - (1-p).^nrI;
        
        % Veletlenszeru megfertozes
        newInfection = rand(size(csucs)) < pInf & isS(csucs);
        
        % Veletlenszeru meggyogyulas
        newRecovery = rand(size(csucs)) < q & isI(csucs);
        
        % Veletlenszeru meggyogyulas
        newImmunityLoss = rand(size(csucs)) < w & isR(csucs);
        
        % Ha mar semmi valtozas sem tortent a populacioban, akkor hagyja abba
        if sum(isI(csucs)) == 0
            break
        end
        
        % Frissitjuk a csucsok erteket
        csucs(newInfection) = 1;
        csucs(newRecovery) = 2;
        csucs(newImmunityLoss) = 0;
        
        fert_num(simNr,i+1)=sum(isI(csucs))/G.nv; %fert�z�ttek ar�nya
        set(h(isS(csucs)),'MarkerFaceColor',[0 1 0]); %cs�csok sz�nez�se
        set(h(isI(csucs)),'MarkerFaceColor',[1 0 0]);
        set(h(isR(csucs)),'MarkerFaceColor',[0 0 1]); %%%%%
        Pl(simNr).YData = fert_num(simNr,:);
        drawnow
    
        % waitforbuttonpress;
        % pause(0.1)
    end

end