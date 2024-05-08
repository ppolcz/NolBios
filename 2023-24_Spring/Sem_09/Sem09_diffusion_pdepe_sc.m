%% Válogatott modellek a matematikai biológiában

%% 2. gyakorlat
% Hővezetés egyenlete végtelen rúdban
% Diffúzió egyenlet megoldása a beépített |pdepe függvény segítségével

%% Hővezetés egyenlete adott kezdeti feltétel mellett
% Hővezetés egyenlete: $u'_t = D u''_{xx}$, ahol $k=1$.
% 
% A Matlab |pdepe| függvénye a következő típusú egyváltozós időfüggő PDE-t tudja 
% megoldani: 
% 
% $$    c(x,t,u,u'_x) \, u'_t = x^{-m} \frac{\partial}{\partial    x} \Big(x^m 
% \, f(x,t,u,u'_x)\Big) + s(x,t,u,u'_x),$$
% 
% $p(x,t,u) + q(x,t) \, f(x,t,u,u'_x) = 0$ típusú peremfeltétellel.
% 
% ahol $m=0,1,2$. 
% 
% Nézzünk néhány példát:
% 
% $$\left\{\begin{array}{l} m = 0 \\ c(x,t,u,u'_x) = 0 \\f(x,t,u,u'_x) = D \, 
% u'_x \\s(x,t,u,u'_x) = 0 \\p(x,t,u) = 0 \\q(x,t) = 1\end{array} \right. \Rightarrow\left\{\begin{array}{ll}u'_t 
% = D u''_{xx} & \text{Diffúzió egyenlete}\\u'_x(0,t) = u'_x(L,t) = 0 & \text{Neumann 
% perem}\end{array}\right.$$
% 
% Vagy:
% 
% $$\left\{\begin{array}{ll} m = 0 \\ c(x,t,u,u'_x) = 0 \\f(x,t,u,u'_x) = D 
% \, u'_x \\s(x,t,u,u'_x) = u - u^2 \\p(x,t,u) = u-a & \text{b.o.}\\p(x,t,u) = 
% u-b & \text{j.o.}\\q(x,t) = 0\end{array} \right. \Rightarrow\left\{\begin{array}{ll}u'_t 
% = D u''_{xx} + cu(1 - u) & \text{Diffúzió egyenlete gerjesztéssel}\\u(0,t) = 
% a & \text{Dirichlet perem} \\u(L,t) = b\end{array}\right.$$
% 
% A |pdepe(m,pde,ic,bc,x,t)| a következő argumentumokat várja:

%% 
% # m = 0, 1, 2 (egyértelmű)
% # pde: $(x,t,u,u'_x) \mapsto \Big( c(x,t,u,u'_x),\, f(x,t,u,u'_x),\, s(x,t,u,u'_x)\Big)$. 
% Tehát egy olyan függvény, aminek négy argumentuma van $(x,t,u,u'_x)$ és három 
% visszatérési értéke: $\Big( c(x,t,u,u'_x),\, f(x,t,u,u'_x),\, s(x,t,u,u'_x) 
% \Big)$. Ez határozza meg a PDE alakját.
% # ic: $x \mapsto g(x)$ kezdeti feltétel
% # bc: $\Big( x_\text{left},u(x_\text{left},t),x_\text{right}, u(x_\text{right},t), 
% t \Big) \mapsto \Big( p_\text{left}, q_\text{left}, p_\text{right}, q_\text{right} 
% \Big)$ peremfeltétel
% # x: felosztás térdimenzió mentén
% # t: felosztás idő dimenzió mentén
% PDE megadása
% A c, f, s függvényeket egy olyan függvény segítségével deklarálom, aminek 
% négy argumentuma van $(x,t,u,u'_x)$, és három értéket ad vissza: $c(x,t,u,u'_x),\, 
% f(x,t,u,u'_x),\, s(x,t,u,u'_x)$, lásd pl. |pde_diffusion_1| függvényt. 
% 
% Ezt a függvényt "inline" is lehet deklarálni egy anonim függvény segítségével 
% (|f = @(x) sin(x)|). Annak érdekében, hogy az anonim függvény több értéket is 
% vissza tudjon adni felhasználom a |deal| függvényt, az úgymond "dealert", amely 
% kiosztja az argumentumában található értékeket különálló visszatérési értékekként.
% Kezdeti feltétel megválasztása
% $$u(x,0) = g(x)$$

g = @(x) sign(sin(pi*x)) + sign(cos(3*pi*x)) + sign(cos(10*pi*x + 0.2));
% g = @(x) sign(sin(2*pi*x));
g = @(x) rand(size(x))*2-1;

x_lim = [0 2.5];
x = linspace(x_lim(1),x_lim(2),1000);

figure, plot(x,g(x),'LineWidth',3)

%% 
% 
% PDE megoldása
% Idő szerinti felosztás megadása:

% t = linspace(0,1,00);
t = logspace(-5,1,300);

%% 
% Végül a megoldás (az x szerinti felosztást már korábban megadtuk):

D = 1;
pde = @(x,t,u,ux) deal( 1 , 0.01*D*ux , 10*u*(1-u^2));

m = 0;
u = pdepe(m,pde,g,@bc_Neumann,x,t);

% Megjelenítés
% Amikor időben több hasonló ábrát szeretnénk egymás után kirajzolni úgy mint 
% egy animációt, akkor csak az első időpillanathoz tartozó ábrákat rajzolom ki 
% (|plot|, |surf|, |mesh|, stb... függvényekkel). Majd az animáció során a már 
% meglévő ábrákat frissítem. Így jóval gyorsabb a kirajzolás, kevésbé fog szaggatni.
% 
% Létrehozok egy Figure ablakot Visibile = 'on' argumentummal, annak érdekében, 
% hogy az ábra kiugorjon a LiveScript-ből. Ennek hiányában az ábra a LiveScript-ben 
% inline fog megjelenni. Én általában módosítani szoktam az ábra méretét is.

fig = figure(1231);
fig.Position(3:4) = [1111 390];

%% 
% Subplot-ot helyett a |tiledlayout|-ot javaslom. Könnyebb használni.

Tl = tiledlayout(2,2,"TileSpacing","compact","Padding","tight");

%% 
% Az első ábrán a teljes $u(x,t)$ függvényt ábrázolom 3D-ben.

nexttile([2,1]);
surf(x,t,u), light, shading interp
xlabel('Distance x')
ylabel('Time t')

%% 
% A második ábrán $v(x) = u(x,t_i)$ függvényt ábrázolom minden egyes $t_i$ időpillanatban. 
% Ábrázolom továbbá $v(x)$ függvény csúcsait és gyökeit. Ebből egy animáció fog 
% kikerekedni.

nexttile; hold on, grid on, box on;
Pl = plot(x,g(x));
ylim([-3,3]) 

%% 
% A kezdeti feltétel függvényére a csúcsokat és gyököket meg sem próbálom kirajzolni, 
% ámde az ábrák objektumait létrehozom és az ezekre mutató pointert eltárolom, 
% hogy az animáció során befolyásolni tudjam.

Pl_pPeaks = plot(0,0,'^');
Pl_nPeaks = plot(0,0,'v');
Pl_Roots = plot(0,0,'k.','MarkerSize',15);
Tit2_Str = @(t) "$v(x) = u(x,t)$, ahol $t = " + sprintf('%9.7f',t) + "$ [sec]";
Tit2 = title(Tit2_Str(0),'Interpreter','latex','FontSize',13);

%% 
% Ábrázolni szeretném a függvény gyökeit is egy külön ábrán.

nexttile; hold on, box on, grid on;
Pl_tRoots = plot(0,NaN,'k.','MarkerSize',7);
xlim(x_lim);
ylim(log10(t([1,end])));
title('Gyökök mozgása')
drawnow

%% 
% Az animációs során minden időpillanatban 

for i = 1:numel(t)
    v = u(i,:);

    Pl.YData = v; 

    % Megkeresem a lokális maximumokat
    [PKS,LOCS] = findpeaks(v,x);
    % Frissítem az adatpontokat
    Pl_pPeaks.XData = LOCS;
    Pl_pPeaks.YData = PKS;

    % Megkeresem a lokális minimumokat
    [PKS,LOCS] = findpeaks(-v,x);
    % Frissítem az adatpontokat
    Pl_nPeaks.XData = LOCS;
    Pl_nPeaks.YData = -PKS;

    % Megkeresem a gyököket
    LOCS = findroots(v,x);
    % Frissítem az adatpontokat
    Pl_Roots.XData = LOCS;
    Pl_Roots.YData = LOCS*0;

    % Hozzáadom a gyökök adatpontjai a gyökök időbeni mozgását leíró ábrához
    Pl_tRoots.XData = [Pl_tRoots.XData LOCS];
    Pl_tRoots.YData = [Pl_tRoots.YData LOCS*0+log10(t(i))];

    % Frissítem a második ábra címét, hogy mindig az aktuális idő jelenleg meg
    Tit2.String = Tit2_Str(t(i));

    drawnow
    pause(0.0001)
end