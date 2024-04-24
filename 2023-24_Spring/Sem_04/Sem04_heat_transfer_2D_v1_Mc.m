%% Heat Equation on a Square Plate
% The basic heat equation is
% 
% $$\frac{\partial u}{\partial t} = k \Delta u$$

geomFileName = "pdeGeom_SQ1-C1_C2.mat";
if ~exist(geomFileName,'file')
    pderect([-1 1 -1 1])
    pdecirc(-0.5,-0.5,0.25)
    pdecirc(0.5,0.5,0.25)
    pdeModeler
    save(geomFileName,"gd","sf","ns")
else
    load(geomFileName,"gd","sf","ns")
end

%%

pdeGeom = decsg(gd,sf,ns);
pdegplot(pdeGeom,"EdgeLabels","on","FaceLabels","on")
axis equal

%% PDE model
% Create a PDE Model with a single dependent variable

numberOfPDE = 1;
pdeModel = createpde(numberOfPDE);
geometryFromEdges(pdeModel,pdeGeom);

%% Problem Definition

k = 0.025;

c = k;
a = 0;
f = 0;
d = 1;

specifyCoefficients(pdeModel,"m",0,...
                             "d",1,...
                             "c",1,...
                             "a",0,...
                             "f",0);
%% Apply Boundary Conditions

% Solution is zero at all four outer edges of the square
% applyBoundaryCondition(pdeModel,'Edge',1:size(pdeGeom,2), 'u', 0);

applyBoundaryCondition(pdeModel,"neumann", ...
    "Edge",1:8, ...
    "g",0, ...
    "q",0);

%% Generate Mesh

msh = generateMesh(pdeModel,'Hmax',0.1);
figure;
pdemesh(pdeModel);
axis equal

%% Initial Conditions

setInitialConditions(pdeModel,0);
setInitialConditions(pdeModel,1,"Face",2);

%%

tlist = logspace(-2,1,150);
results = solvepde(pdeModel,tlist);

sol = results.NodalSolution;

%%

fig = figure(123);
fig.Visible = 'on';
delete(fig.Children);
ax = axes(fig);

Tri = msh.Elements(1:3,:)';
Sf = trisurf(Tri,msh.Nodes(1,:)',msh.Nodes(2,:)',sol(:,1));

for i = 1:numel(tlist)
    Sf = trisurf(Tri,msh.Nodes(1,:)',msh.Nodes(2,:)',sol(:,i));
    zlim([0,1])
    clim([0,1])
    drawnow
end

%%
return
%% Idointervallum, lepeskoz

nframes = 50;
tlist = linspace(0,1,nframes);
%% Itt tortenik a PDE megoldasa

u1 = parabolic(u0,tlist,pdeModel,c,a,f,d);
%% Plot FEM Solution
% To speed up the plotting, we interpolate to a rectangular grid.

fig = figure('Position', [168 562 982 386], 'Color', 'white');
colormap(hot);
x = linspace(-1,1,301);
y = x;
[~,tn,a2,a3] = tri2grid(p,t,u0,x,y);
umax = max(max(u1));
umin = min(min(u1));

for j = 1:nframes,
    u = tri2grid(p,t,u1(:,j),tn,a2,a3);

    figure(fig)
    subplot(121);
    surf(x,y,u);
    view([0,89.999])
    caxis([umin umax]);
    axis([-1 1 -1 1 0 2]);
    shading interp;

    subplot(122);
    surf(x,y,u);
    view([10,60])
    caxis([umin umax]);
    axis([-1 1 -1 1 0 1.2]);
    shading interp
    light
    title(sprintf('time = %g', tlist(j)))

    pause(0.05)

    if j == 20
        drawnow
        persist.savefig(fig, 'heat_diffusion_poster.png');
        snapnow
        pause(0.5)
    end
    
    % Felvétel
    % frames(j) = getframe(fig);
end

% v = VideoWriter(persist.timefl('fig','heat_diffusion.avi'));
% open(v)
% writeVideo(v,frames)
% close(v)