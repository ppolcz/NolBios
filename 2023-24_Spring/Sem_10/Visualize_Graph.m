function [Pl,Tl,ax,fig] = Visualize_Graph(G,FigNr)
arguments
    G,
    FigNr = 1
end
%%
%  Author: Peter Polcz (ppolcz@gmail.com) 
%  Created on 2024. May 14. (2023a)
%

fig = figure(FigNr);
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

end
