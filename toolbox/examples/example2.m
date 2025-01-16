%% Create DEM
% Here we use the Grindelwald to get a grid geometry, but we built an
% inclined plane that abruptly transitions into a flat area
DEM = GRIDobj("data\grindelwald.txt");
DEM = GRIDobj(DEM);
DEM.Z(:,600:end) = repmat(linspace(0,5000,DEM.size(2)-600+1),DEM.size(1),1);

% The detached height
H   = GRIDobj(DEM);
H.Z(300:400,900:1000) = 100;

[H2,T] = ttminvoellmy(DEM,H,"video",true,'maxtime',500,'camorbit',[0 0], ...
    'controls',false,'mu',0.75,'xi',500,'videoname','vid.mp4','hmin',0.001,...
    'vc',4);
