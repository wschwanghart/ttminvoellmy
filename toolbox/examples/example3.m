%% Create DEM
% Here we use the Grindelwald to get a grid geometry, but we built an
% inclined plane that abruptly transitions into a flat area
DEM = GRIDobj("data\grindelwald.txt");
DEM = GRIDobj(DEM);
DEM = DEM + repmat(linspace(0,300,DEM.size(2)),DEM.size(1),1);
DEM = DEM + repmat(linspace(-1,1,DEM.size(1))',1,DEM.size(2)).^2 * 2500;

% The detached height
H   = GRIDobj(DEM);
H.Z(100:200,800:900) = 150;

[H2,T] = ttminvoellmy(DEM,H,"video",true,'maxtime',500,'camorbit',[180+135 -20], ...
    'controls',false,'mu',0.2,'xi',200,'videoname','vid.mp4','hmin',0.01,...
    'vc',4,'clim',[0 100],'maxdt',0.05);
