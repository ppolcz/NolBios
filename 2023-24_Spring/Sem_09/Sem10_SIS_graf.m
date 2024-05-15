
rajz = 1;
fokszam = 1;
n = 1000;
n0 = 10;
L = 2;
L0 = 5;
[G] = BAmodel_v2(n,n0,L,L0);

fig = figure(1);
Tl = tiledlayout(1,1,"Padding","none");
ax = nexttile;
Pl = plot(G,'Layout','force3','Iterations',10,'UseGravity',true);
% Pl = plot(G,'Layout','subspace3','Dimension',10);
Pl.NodeCData = log10(degree(G)+1);
Pl.MarkerSize = 5;
set(ax,"Visible","off","Clipping","off")
set(fig,"Color",[1,1,1])
colormap parula
rotate3d on
