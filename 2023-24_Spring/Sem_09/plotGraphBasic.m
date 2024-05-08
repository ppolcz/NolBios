function plotGraphBasic(G,markerSize,addText)
% gráfrajzoló rutin, kör alakban felrajzolja a gráfot.
% Created by Pablo Blinder. blinderp@bgu.ac.il
% Updated by David Tekan. tekda@digitus.itk.ppke.hu

figure;

% Gráfrajzoló subrutin. Elhelyezi a gráf nv csúcsát egy körvonal mentén
center=[0,0];
theta=linspace(0,2*pi,G.nv+1);
rho=ones(1,G.nv+1);
[X,Y] = pol2cart(theta',rho');
X=X+center(1);
Y=Y+center(2);
x=(X(1:end-1)*10)';
y=(Y(1:end-1)*10)';
%

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
kv(kv<1)=1;
Pv=num2cell(map(kv,:),2);
if kvGroups==1; kvGroups=2; end 
set(gca,'Clim',[1 max(kvGroups)]);


Pn(1)={'MarkerFaceColor'};


h = nan(G.nv,1);
for i=1:G.nv
    h(i,:)=plot(x(i),y(i),'ko');
end


if addText
    for i=1:G.nv
       text(x(i)+0.1*x(i),y(i)+0.1*y(i),num2str(i));
    end
end

set(h,'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerSize',markerSize,Pn,Pv);
set(gca,'Visible','Off','YDir','reverse');
colormap(map);
hc=colorbar;

set(hc,'FontSize',8,'FontW','Demi')
set(hc,'Visible','off')
set(gcf,'Color','w')