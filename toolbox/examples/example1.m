%% Run Stefan Hergarten's minvoellmy model


%% Load data
DEM = GRIDobj("data\grindelwald.txt");

%% Create detached volume
pos = [...
    6.4348    1.5869
    6.4360    1.5879
    6.4370    1.5885
    6.4387    1.5891
    6.4395    1.5897
    6.4392    1.5912
    6.4371    1.5933
    6.4359    1.5928
    6.4342    1.5908
    6.4319    1.5891] * 1e5;

MS = struct('Geometry','Polygon',...
    'X',pos(:,1),'Y',pos(:,2));

LS = polygon2GRIDobj(DEM,MS,'waitbar',false);
H  = LS*20;

%% Run minvoellmy

[H2,T] = ttminvoellmy(DEM,H,"gif",false,'maxtime',200,'camorbit',[160 0]);

