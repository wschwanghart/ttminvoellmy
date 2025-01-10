function [H,T] = ttminvoellmy(DEM,H,options)
%TTMINVOELLMY Rock avalanche model 
%
% Syntax


arguments
    DEM   % GRIDobj of DEM: All units must be in [m]
    H     % GRIDobj of thickness of failed mass [m]
    options.mu = 0.2
    options.xi = 500
    options.vc (1,1) {mustBeNumeric,mustBeNonnegative} = 4   % crossover velocity [m s^-1]
    options.hmin = 0.01
    options.dmin = 0
    options.cent = true
    options.g  (1,1) {mustBeNumeric}  = 9.81
    options.maxsteps = inf
    options.maxtime  = 1000 % [s]
    options.maxdt    = 1% [s]
    options.cfl      = 0.7

    options.plot     (1,1) = true
    options.plottype = 'surf'
    options.plotmask = ~GRIDobj(DEM,'logical')
    options.camorbit (1,2) {mustBeNumeric} = [0,0] 
    options.gif      (1,1) = false
    options.gifname  = 'minvoellmy.gif'
    options.colormap = flipud(ttscm('batlowW',255,[20 100]))
    options.clim     (1,2) {mustBeNumeric} = [0, 5]
    options.controls = true
    options.output (1,1) {mustBeNumeric,mustBeNonnegative} = inf
    options.position (1,4) = [0.20 .1 .70 .80]

end

% Do some input checking
validatealignment(DEM,H)
b = DEM.Z;
h = H.Z;

% Extract cell sizes
dx = DEM.cellsize;
dy = DEM.cellsize;

% Initiate the simulation class
o = minvoellmy2d(b,h,dx,dy,...
    options.mu,options.xi,options.vc,...
    options.hmin,options.dmin,options.cent,options.g);

% Preallocate array for flow depth
H  = DEM;

% counter
counter = 1;
% track time
t = 0;

% DEMp for plotting
DEMp = crop(DEM,options.plotmask);

% Set up figure and controls
if options.plot
    if options.controls
        fhcontrol = uifigure('Units','Normalized',...
            'Position',[0 .60 .20 .30],...
            'Name','Controls',...
            'NumberTitle','off','WindowStyle','alwaysontop');

        g = uigridlayout(fhcontrol,[4 1]);
        g.RowHeight = {'fit','1x','1x','1x'};
        g.ColumnWidth = {'1x'};
        htimetext   = uitextarea(g,"Value",num2str(0));
        hpausebutton = uibutton(g,"Text","Pause simulation",...
            "ButtonPushedFcn",@(btn,event) pauseOrContinueLoop());
        hstopbutton = uibutton(g,"Text","Stop simulation",...
            "ButtonPushedFcn",@(btn,event) stopLoop());
       
        

        % Create a variable to control the loop
        running = true;
        ispaused = false;
        drawnow

        % Create another axis with the figure
        fh  = figure('Units','Normalized','Position',options.position);
        if options.gif
            % If a gif will be produced, then the figure should not be
            % resizable.
            fh.Resize = false;
        end
        ax  = axes("Parent",fh);
        ax.Position = [0 0 1 1];
        H.Z = o.h;
        Hp   = crop(H, options.plotmask);
        h    = surf(DEMp + Hp,Hp,'block',true,'parent',ax);
        camorbit(ax,options.camorbit(1),options.camorbit(2))
        colormap(ax,options.colormap)
        camlight

        % colorbar
        clim(options.clim)
        lastcamorbit = options.camorbit;

        axis off
        material dull

        title(gca,num2str(t))
        drawnow;

        % Gif
        if options.gif
            gif(options.gifname, 'frame', fh,'overwrite',true)
        end

        hslider = uislider(g,'Value',options.camorbit(1),'Limits',[0 360],...
            'ValueChangingFcn',@(src,event,xr) setCamorbit());
    end

else
    running = true;
    ispaused = false;
end

T = table(0,H,'VariableNames',{'Time','H'});

% Run loop    
while counter <= options.maxsteps && t <= options.maxtime && running

    dt = step(o,options.maxdt,options.cfl);
    t  = t+dt;
    
    if options.controls
    htimetext.Value = num2str(round(t,1));
    end
    % V.Z = sqrt(o.uh.^2 + o.vh.^2);

    while ispaused
        pause(0.2)
    end

        
    if counter == 1 || mod(counter,10) == 0
        H.Z = o.h;
        Hp   = crop(H, options.plotmask);

        h.ZData = DEMp.Z + Hp.Z;
        h.CData = Hp.Z;

        title(gca,num2str(t))
        drawnow; 
    end

    if options.gif && mod(counter,10) == 0
        if options.gif
            gif
        end
    end

    if ~isinf(options.output) && ((mod(counter,options.output)) == 0)
        H.Z = o.h;
        T   = [T;table(t,H,'VariableNames',{'Time','H'})];
    end

    % u = max(o.uh(:));
    counter = counter +1 ;
end

H.Z = o.h;

if options.controls
    close(fhcontrol)
end
% if options.gif
%     gif('clear')
% end


% Nested functions 
    function stopLoop()
        running = false;
    end
    function pauseOrContinueLoop()
        if ~ispaused 
            ispaused = true;
            hpausebutton.Text = 'Resume simulation';
        else
            ispaused = false;
            hpausebutton.Text = 'Pause simulation';
        end
    end
    function setCamorbit()
        co = lastcamorbit;
        val = hslider.Value;
        newvalue = val-co(1);
        camorbit(ax,newvalue,0)
        lastcamorbit(1) = newvalue;
    end

end
