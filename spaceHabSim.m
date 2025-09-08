function spinStation_minimal
% --- pre-run scrub: kill any leftover figures/timers from prior runs ---
set(0,'ShowHiddenHandles','on');
try
    delete(findall(0,'Type','figure','Tag','StationSim'));
    oldTimers = timerfindall('Tag','StationSimTimer');
    if ~isempty(oldTimers), stop(oldTimers); delete(oldTimers); end
catch
end
set(0,'ShowHiddenHandles','off');

%% ---- knobs ----
R_ring      = 20;       % OUTER walking radius [m]
nSpokes     = 6;        % spokes
nU          = 64;       % torus resolution (along)
nV          = 48;       % torus resolution (around)
spokeFacets = 28;       % cylinder facets
nStars      = 2500;      % star count

% Average walking speed used for Coriolis readout
WALK_SPEED = 1.40;  % m/s

% --- rotation input in RPM (sign sets direction) ---
spinRPM     = 7;                        % rev/min
spinRateSec = spinRPM * (2*pi/60);      % rad/s (kept for dependencies)
period      = 0.03;                     % ~33 Hz timer

% ---- UI theme (shared with bottom-right controls) ----
MENU_BG  = [0.05 0.06 0.08];
MENU_FG  = [0.85 0.90 1.00];
MENU_ACC = [0.30 0.60 1.00];
MENU_FONT= 'Consolas';

%% ---- station size policy (R_ring is OUTER radius; tube Ø fixed at 2 m) ----
R_OUTER_MIN = 2.0;          % m
R_OUTER_MAX = 3.0 * 60.0;   % m  (≈ 180 m, ~3× ISS scale)
R_outer     = max(R_OUTER_MIN, min(R_ring, R_OUTER_MAX));
R_tube      = 1.0;                          % 2 m diameter => 1 m radius
R_ring_center = R_outer - R_tube;           % torus centerline radius

%% ---- figure / axes (two-layer: background & foreground) ----
opengl hardware;
fig = figure('Name','Rotating Station (HQ+Smooth)','Color','k', ...
             'NumberTitle','off','Renderer','opengl','Interruptible','off', ...
             'Tag','StationSim');
try fig.GraphicsSmoothing = 'on'; end

% BACKGROUND axes: planet + stars
axBG = axes('Parent',fig,'SortMethod','depth','Color','k');
axis(axBG,'equal'); axis(axBG,'off'); axis(axBG,'vis3d');
hold(axBG,'on'); view(axBG,35,20); camproj(axBG,'perspective');
camtarget(axBG,[0 0 0]); camva(axBG,12);

% FOREGROUND axes: station
axFG = axes('Parent',fig,'SortMethod','depth','Color','none'); % transparent
axis(axFG,'equal'); axis(axFG,'off'); axis(axFG,'vis3d');
hold(axFG,'on'); view(axFG,35,20); camproj(axFG,'perspective');
camtarget(axFG,[0 0 0]); camva(axFG,12);

% Link cameras so motion is unified
lp = linkprop([axBG, axFG], ...
    {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle','Projection'});
setappdata(fig,'CamLink',lp);
uistack(axFG,'top');

try disableDefaultInteractivity(axBG); end
try disableDefaultInteractivity(axFG); end
set([axBG axFG fig],'HitTest','off');

%% ---- station group (on axFG) ----
G = hgtransform('Parent',axFG);

% Torus (ring) using centerline radius
[U,Vv] = meshgrid(linspace(0,2*pi,nU), linspace(0,2*pi,nV));
[Xt,Yt,Zt] = torusMesh(R_ring_center, R_tube, U, Vv);
torusSurf = surf(axFG, Xt, Yt, Zt, ...
    'EdgeColor','none', 'FaceLighting','gouraud','FaceAlpha',1.0, ...
    'AmbientStrength',0.25,'DiffuseStrength',0.7, ...
    'SpecularStrength',0.5,'SpecularExponent',25, ...
    'HitTest','off','PickableParts','none', ...
    'FaceColor',[0.8 0.8 0.85], 'Parent', G);

% Spokes (to centerline) — fixed 0.5 m radius
spokeHs = gobjects(0);
if nSpokes>0
    th = linspace(0,2*pi,nSpokes+1); th(end) = [];
    spokeRad = 0.9;   % meters (fixed)
    for i = 1:nSpokes
        p0 = [0 0 0];
        p1 = [R_ring_center*cos(th(i)) R_ring_center*sin(th(i)) 0];
        [SX,SY,SZ] = cylBetween(p0,p1,spokeRad,spokeFacets);
        s = surf(axFG, SX,SY,SZ, 'EdgeColor','none', ...
            'FaceLighting','gouraud', ...
            'AmbientStrength',0.25,'DiffuseStrength',0.7, ...
            'SpecularStrength',0.5,'SpecularExponent',25, ...
            'HitTest','off','PickableParts','none', ...
            'FaceColor',[0.8 0.8 0.85]);
        set(s,'Parent',G);
        spokeHs(end+1,1) = s; %#ok<AGROW>
    end
end

% Hex hub
silver = [0.8 0.8 0.85];
hubHexRadius = 0.25 * R_ring_center;
hubDepth     = 2.2 * R_tube;
[Fhub, Vhub] = regularPrism(6, hubHexRadius, hubDepth);
hubH = patch('Parent', G, 'Faces', Fhub, 'Vertices', Vhub, ...
    'EdgeColor','none', 'FaceColor', silver, ...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.25, 'DiffuseStrength',0.7, ...
    'SpecularStrength',0.5, 'SpecularExponent',25, ...
    'HitTest','off','PickableParts','none');

% Foreground limits (around station)
Rmax = 200;
axis(axFG,[-Rmax Rmax -Rmax Rmax -0.6*Rmax 0.6*Rmax]);

% Station light (infinite dir; station only)
sunDir = [-1 -0.2 0.3]; sunDir = sunDir / norm(sunDir);
light(axFG,'Position',sunDir,'Style','infinite','HitTest','off');

%% ---- planet (on axBG, live-lit) ----
SCALE_FAR   = 1/100;     % meters → “km units” for far objects
Rplanet_m   = 6371000;    % Earth radius (m)
orb_alt_m   = 200000;     % 200 km
planetCtr_m = [Rplanet_m + (-20*orb_alt_m), ...
               Rplanet_m + ( 10*orb_alt_m), ...
              -(Rplanet_m + ( -3*orb_alt_m))];

Rplanet   = Rplanet_m   * SCALE_FAR;
planetCtr = planetCtr_m * SCALE_FAR;

baseTex = proceduralEarthTexture(1024, 2048);

[xe,ye,ze] = sphere(64);
Xe = Rplanet*xe + planetCtr(1);
Ye = Rplanet*ye + planetCtr(2);
Ze = Rplanet*ze + planetCtr(3);

planet = surf(axBG, Xe, Ye, Ze, ...
    'EdgeColor','none','FaceColor','texturemap', ...
    'CData', baseTex, 'CDataMapping','direct', ...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.10,'DiffuseStrength',0.78, ...
    'SpecularStrength',0.20,'SpecularExponent',4, ...
    'BackFaceLighting','unlit', ...
    'HitTest','off','PickableParts','none');
set(planet,'Clipping','off');

% Planet/sun light on background
oldSun = findall(axBG,'Type','light','Tag','StationSun');
if ~isempty(oldSun), delete(oldSun); end
light(axBG,'Style','infinite','Position',sunDir(:).', ...
          'Color',[1 1 1],'Tag','StationSun');

%% ---- starfield (on axBG, large fixed shell) ----
rng(45);
dirs = randn(nStars,3); dirs = dirs ./ vecnorm(dirs,2,2);
radJitter = 0.9 + 0.2*rand(nStars,1);
Rstar = 2 * (norm(planetCtr) + Rplanet);
S = (dirs .* radJitter) * Rstar;

starsH = scatter3(axBG, S(:,1), S(:,2), S(:,3), 6, 'filled', ...
                  'MarkerFaceColor',[0.85 0.9 1], 'MarkerEdgeColor','none', ...
                  'HitTest','off','PickableParts','none');
set(starsH,'Clipping','off');

% Sun sprites (on axBG, once)
Rsun   = 2 * (norm(planetCtr) + Rplanet);
sunPos = Rsun * sunDir(:).';
sizes  = (0.0001*(norm(planetCtr) + Rplanet))*[1400 900 600 380];
cols   = [0.35 0.35 0.40;
          0.60 0.65 0.75;
          0.85 0.90 1.00;
          1.00 1.00 1.00];
for k = 1:numel(sizes)
    h = scatter3(axBG, sunPos(1), sunPos(2), sunPos(3), sizes(k), ...
                 'o','filled','MarkerFaceColor',cols(k,:), ...
                 'MarkerEdgeColor','none','HitTest','off','PickableParts','none');
    set(h,'Clipping','off');
end

% Background axes limits to include planet + stars
Rscene = norm(planetCtr) + Rplanet + Rstar*1.1;
axis(axBG, [-Rscene Rscene -Rscene Rscene -Rscene Rscene]);

%% ---- HUDs ----
fpsBox = annotation(fig,'textbox',[0.01 0.95 0.2 0.04], ...
    'String','FPS: --','FontSize',10,'Color',[0.9 0.95 1], ...
    'FontName','Consolas','EdgeColor','none','BackgroundColor','none', ...
    'HorizontalAlignment','left','VerticalAlignment','top');

% Taller box to fit extra lines (Δg + Coriolis)
gHUD = annotation(fig,'textbox',[0.01 0.01 0.40 0.13], ...
    'String', sprintf(['g_{head}: -- g\n' ...
                       'g_{feet}: -- g\n' ...
                       '\\Delta g: -- g\n' ...
                       'Coriolis @ %.2f m/s (feet): with: -- g   against: -- g'], ...
                       WALK_SPEED), ...
    'FontSize',10,'Color',[0.9 0.95 1], 'FontName','Consolas', ...
    'EdgeColor','none','BackgroundColor','none', ...
    'HorizontalAlignment','left','VerticalAlignment','bottom', ...
    'Interpreter','tex');

%% ---- TOP-RIGHT DROPDOWN MENU ----
% Hamburger (≡) button
menuBtn = uicontrol('Parent',fig,'Style','pushbutton','Units','normalized', ...
    'Position',[0.955 0.948 0.035 0.045], 'String', char(8801), ... % '≡'
    'FontName', MENU_FONT, 'FontSize', 16, ...
    'ForegroundColor', MENU_FG, 'BackgroundColor', [0.12 0.13 0.16], ...
    'TooltipString','Menu','Callback',@toggleMenu);

% Dropdown container (hidden by default). Styled like control panel.
menuPanel = uipanel('Parent',fig,'Units','normalized', ...
    'Position',[0.74 0.68 0.25 0.28], 'Visible','off', ...
    'BackgroundColor', MENU_BG, 'BorderType','line', ...
    'ForegroundColor', MENU_ACC, 'Title','  Menu  ', ...
    'FontName', MENU_FONT, 'FontSize',10,'HighlightColor',MENU_ACC);

% Child panels for each level (stacked; only one visible at a time)
menuTop     = uipanel('Parent',menuPanel,'Units','normalized','Position',[0 0 1 1], ...
                 'BackgroundColor',MENU_BG,'BorderType','none','Visible','on');
menuAbout   = uipanel('Parent',menuPanel,'Units','normalized','Position',[0 0 1 1], ...
                 'BackgroundColor',MENU_BG,'BorderType','none','Visible','off');
menuOutputs = uipanel('Parent',menuPanel,'Units','normalized','Position',[0 0 1 1], ...
                 'BackgroundColor',MENU_BG,'BorderType','none','Visible','off');

% Reusable button factory
    function b = mkBtn(parent, pos, label, cb)
        b = uicontrol(parent,'Style','pushbutton','Units','normalized', ...
            'Position',pos,'String',label,'FontName',MENU_FONT, ...
            'ForegroundColor',MENU_FG,'BackgroundColor',[0.15 0.16 0.20], ...
            'Callback',cb);
    end

% Reusable small "Back" button
    function mkBack(parent, cb)
        mkBtn(parent,[0.04 0.86 0.24 0.10],'◀ Back',cb);
    end

% -------- Level 0: Top --------
% Layout: tall buttons stacked
mkBtn(menuTop,[0.04 0.68 0.92 0.18],'About',      @(~,~) showPanel(menuAbout));
mkBtn(menuTop,[0.04 0.46 0.92 0.18],'Inputs',     @(~,~) showText('Inputs', ...
    sprintf('Station radius and RPM.\n\nUse the lower-right UI to change them.')));
mkBtn(menuTop,[0.04 0.24 0.92 0.18],'Outputs',    @(~,~) showPanel(menuOutputs));
mkBtn(menuTop,[0.04 0.02 0.92 0.18],'Controls',   @(~,~) showText('Controls', ...
    sprintf(['Drag-click to swivel (orbit camera).\n' ...
             'Scroll wheel to zoom.\n' ...
             'Use the lower-right UI to control station radius and RPM.'])));

% -------- Level 1: About --------
mkBack(menuAbout, @(~,~) showPanel(menuTop));
mkBtn(menuAbout,[0.04 0.68 0.92 0.18],'Simulation design', @(~,~) showText('Simulation design', ...
    sprintf(['This is a simple rotating space habitat simulator designed to help\n' ...
             'gauge station size and RPM constraints for human comfort and engineering limits.\n\n' ...
             'The outer ring is set to 2 meters in diameter (tube Ø). Imagine a person\n' ...
             'standing with feet on the outer diameter and head at the inner diameter;\n' ...
             'this gives a clear visual reference while you vary radius and RPM.\n\n' ...
             'Environment note: background stars, planet, and sun are purely ambiance for scale.'])));
mkBtn(menuAbout,[0.04 0.46 0.92 0.18],'Environment', @(~,~) showText('Environment', ...
    sprintf(['The background stars, planet, and sun are not part of the physics.\n' ...
             'They provide ambiance, a sense of scale, and, frankly, a little coding flex.'])));
mkBtn(menuAbout,[0.04 0.24 0.92 0.18],'Constraints', @(~,~) showText('Constraints', ...
    sprintf(['We limit the radius to 2 m minimum and 180 m maximum to span near-future designs.\n' ...
             'Smaller than 2 m is not human-rated (no one could stand), while the upper bound\n' ...
             'is ~2–3× ISS scale to reflect today''s structural/cost limits.'])));

% -------- Level 1: Outputs --------
mkBack(menuOutputs, @(~,~) showPanel(menuTop));
mkBtn(menuOutputs,[0.04 0.68 0.92 0.18],'Centripetal acceleration (“g-force”)', @(~,~) showText('Centripetal acceleration', ...
    sprintf(['We report “g-force” at the top (head) and base (feet) of a 2 m tall person.\n' ...
             'This is the centripetal acceleration from rotation at those radii.'])));
mkBtn(menuOutputs,[0.04 0.46 0.92 0.18],'Delta-g (feet − head)', @(~,~) showText('Delta-g', ...
    sprintf(['Δg = g_feet − g_head. Use this to keep differences small (e.g., 0.05–0.10 g)\n' ...
             'for comfort and physiological tolerance.'])));
mkBtn(menuOutputs,[0.04 0.24 0.92 0.18],'Coriolis (walking with/against spin)', @(~,~) showText('Coriolis effects', ...
    sprintf(['We estimate the change in felt g at the feet for a person walking tangentially\n' ...
             'with or against station spin at a nominal speed. Lower is better for comfort.'])));

% -------- Helpers for menu behavior --------
    function toggleMenu(~,~)
        if strcmp(get(menuPanel,'Visible'),'on')
            set(menuPanel,'Visible','off');
        else
            showPanel(menuTop);
            % anchor dropdown just below the hamburger
            set(menuPanel,'Position',[0.70 0.68 0.28 0.28],'Visible','on');
        end
    end

    function showPanel(p)
        % Hide all level panels, then show target
        set([menuTop menuAbout menuOutputs],'Visible','off');
        set(p,'Visible','on');
        if strcmp(get(menuPanel,'Visible'),'off'), set(menuPanel,'Visible','on'); end
    end

% Centered info overlay (dismiss on any click)
    function showText(titleStr, bodyStr)
        % Hide menu when showing overlay
        set(menuPanel,'Visible','off');
        % Remove previous overlay if any
        old = getappdata(fig,'InfoOverlay');
        if ~isempty(old) && isgraphics(old), delete(old); end
        % Create centered info box
        info = annotation(fig,'textbox','Units','normalized', ...
            'Position',[0.64 0.64 0.35 0.3], ...
            'String', sprintf('%s\n\n%s', titleStr, bodyStr), ...
            'FontName', MENU_FONT, 'FontSize', 12, ...
            'Color', MENU_FG, 'EdgeColor', MENU_ACC, ...
            'BackgroundColor', MENU_BG, ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Interpreter','none');
        setappdata(fig,'InfoOverlay',info);
    end


%% ---- start overlay (centered, disappears on first click) ----
startOverlay = annotation(fig,'textbox', ...
    'Units','normalized', ...
    'Position',[0.15 0.44 0.70 0.12], ...
    'String', sprintf('ROTATING SPACE HABITAT SIMULATOR\nclick & drag anywhere to start'), ...
    'FontSize',16, 'FontName','Consolas', ...
    'Color',[0.95 0.97 1], ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'EdgeColor','none', ...
    'BackgroundColor','none', ...
    'Interpreter','none', ...
    'HitTest','off');
setappdata(fig,'StartOverlay',startOverlay);

%% ---- camera knobs ----
INIT_AZ   = 30;
INIT_EL   = 35;
INIT_CAMVA = R_outer * 0.09;
ZOOM_MIN   = R_outer * 0.006;
ZOOM_MAX   = Rmax   * 0.2;
CAM_UP = [0 0 1];

DRAG_X_SIGN = 1;   % +1 normal, -1 inverted
DRAG_Y_SIGN = -1;   % +1 normal, -1 inverted

view(axFG, INIT_AZ, INIT_EL);
camtarget(axFG, [0 0 0]);
camup(axFG, CAM_UP);
camproj(axFG, 'perspective');
camva(axFG, INIT_CAMVA);
set(axFG,'CameraViewAngleMode','manual');

%% ---- custom camera (orbit on drag; wheel zoom) ----
dragging=false; lastPos=[]; az=INIT_AZ; el=INIT_EL;
fig.WindowButtonDownFcn = @onDown;
fig.WindowButtonUpFcn   = @onUp;
fig.WindowScrollWheelFcn= @onScroll;
updateView();

function onDown(~,~)
    % 1) Dismiss info overlay if present
    infoH = getappdata(fig,'InfoOverlay');
    if ~isempty(infoH) && isgraphics(infoH)
        delete(infoH);
        setappdata(fig,'InfoOverlay',[]);
        return;  % don’t start dragging on this click
    end

    % 2) Close dropdown if open and click is outside it (and not the button)
    if strcmp(get(menuPanel,'Visible'),'on')
        ptFig = get(fig,'CurrentPoint');                % pixels
        panelPix = getpixelposition(menuPanel,true);    % [x y w h]
        btnPix   = getpixelposition(menuBtn,true);
        if ~pointInRect(ptFig, panelPix) && ~pointInRect(ptFig, btnPix)
            set(menuPanel,'Visible','off');
            % don’t return; allow click to start dragging as normal
        end
    end

    % (existing start overlay dismiss)
    h = getappdata(fig,'StartOverlay');
    if ~isempty(h) && isgraphics(h)
        delete(h); setappdata(fig,'StartOverlay',[]);
    end

    if strcmp(get(fig,'SelectionType'),'normal')
        dragging  = true;
        lastPos   = get(fig,'CurrentPoint');
        fig.WindowButtonMotionFcn = @onDrag;
        set(fig,'Pointer','fleur');
    end
end

function onUp(~,~)
    dragging = false;
    fig.WindowButtonMotionFcn = '';
    set(fig,'Pointer','arrow');
end

function onDrag(~,~)
    if ~dragging, return; end
    pos = get(fig,'CurrentPoint'); dp = pos - lastPos; lastPos = pos;
    fpos = getpixelposition(fig);
    sens = 0.1 * camva(axFG);
    az = az - DRAG_X_SIGN * sens * 360 * (dp(1)/max(fpos(3),1));
    el = max(-89, min(89, el + DRAG_Y_SIGN * sens * 180 * (dp(2)/max(fpos(4),1))));
    updateView();
end

function onScroll(~,evt)
    va = camva(axFG);
    va = va * (1 + (evt.VerticalScrollCount)*0.05);
    camva(axFG, max(ZOOM_MIN, min(ZOOM_MAX, va)));
end

function updateView()
    view(axFG, az, el);
    camtarget(axFG,[0 0 0]);
    camup(axFG, CAM_UP);
end

%% ---- CONTROL PANEL (lower-right) ----
ctrl = uipanel('Parent',fig,'Units','normalized','Position',[0.64 0.01 0.35 0.25], ...
    'BackgroundColor',[0.05 0.06 0.08],'BorderType','line','ForegroundColor',[0.3 0.6 1], ...
    'Title',' Controls ','FontName','Consolas','FontSize',10,'HighlightColor',[0.3 0.6 1]);

% Ring Radius
uicontrol(ctrl,'Style','text','Units','normalized','Position',[0.04 0.70 0.45 0.24], ...
    'String','Ring Radius (m)','BackgroundColor',get(ctrl,'BackgroundColor'), ...
    'ForegroundColor',[0.85 0.9 1],'HorizontalAlignment','left','FontName','Consolas');
sRing = uicontrol(ctrl,'Style','slider','Units','normalized','Position',[0.04 0.58 0.72 0.14], ...
    'Min',R_OUTER_MIN,'Max',R_OUTER_MAX,'Value',R_outer,'SliderStep',[0.005 0.05], ...
    'Callback',@onRingSlider);
eRing = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.78 0.56 0.18 0.18], ...
    'String',sprintf('%.2f',R_outer),'BackgroundColor',[0.15 0.16 0.2],'ForegroundColor',[0.95 0.97 1], ...
    'FontName','Consolas','Callback',@onRingEdit);

% RPM (moved up to create extra gap above the buttons)
uicontrol(ctrl,'Style','text','Units','normalized','Position',[0.04 0.36 0.30 0.22], ...
    'String','RPM','BackgroundColor',get(ctrl,'BackgroundColor'), ...
    'ForegroundColor',[0.85 0.9 1],'HorizontalAlignment','left','FontName','Consolas');
sRPM = uicontrol(ctrl,'Style','slider','Units','normalized','Position',[0.04 0.24 0.72 0.14], ...
    'Min',0,'Max',21.15,'Value',spinRPM,'SliderStep',[0.01 0.1], ...
    'Callback',@onRpmSlider);
eRPM = uicontrol(ctrl,'Style','edit','Units','normalized','Position',[0.78 0.22 0.18 0.18], ...
    'String',sprintf('%.2f',spinRPM),'BackgroundColor',[0.15 0.16 0.2],'ForegroundColor',[0.95 0.97 1], ...
    'FontName','Consolas','Callback',@onRpmEdit);

% --- Auto-1g buttons (stay near bottom; now with more space above) ---
btnRPM1g = uicontrol(ctrl,'Style','pushbutton','Units','normalized', ...
    'Position',[0.04 0.06 0.36 0.10], 'String','Set RPM for 1g', ...
    'FontName','Consolas','Callback',@onSetRpm1g);

btnR1g = uicontrol(ctrl,'Style','pushbutton','Units','normalized', ...
    'Position',[0.46 0.06 0.36 0.10], 'String','Set Radius for 1g', ...
    'FontName','Consolas','Callback',@onSetRadius1g);

%% ---- spin + FPS timer (drop late frames, time-based) ----
spinAngle = 0;
t0 = tic; frameTimes = []; winSec = 1.0;
t = timer('ExecutionMode','fixedRate','Period',period, ...
          'BusyMode','drop','TimerFcn',@tick,'StartDelay',0, ...
          'Tag','StationSimTimer');

guard = onCleanup(@() localCleanup(fig, t));
start(t);
fig.CloseRequestFcn = @onClose;

    function tick(~,~)
        nowSec = toc(t0);
        persistent lastSec
        if isempty(lastSec), lastSec = nowSec; end
        dt = nowSec - lastSec; lastSec = nowSec;

        % time-based station spin
        spinAngle = spinAngle + spinRateSec * dt;
        set(G,'Matrix', makehgtform('zrotate', spinAngle));

        drawnow nocallbacks;

        % FPS
        frameTimes(end+1) = nowSec;
        while frameTimes(1) < nowSec - winSec, frameTimes(1) = []; end
        if numel(frameTimes) >= 2
            fps = (numel(frameTimes)-1) / (frameTimes(end)-frameTimes(1));
            set(fpsBox,'String',sprintf('FPS: %4.1f', fps));
        end

        % g at feet/head (2.0 m-tall person on outer rim) + Δg + Coriolis (feet)
        [gFeet, gHead] = stationGs(R_outer, spinRateSec, 2.0);
        deltaG = gFeet - gHead;
        
        % Coriolis magnitude (radial) for tangential walking speed v: |a_cori| = 2*ω*v
        % Convert to g-units using g0 = 9.80665 m/s^2
        gCori = (2 * spinRateSec * WALK_SPEED) / 9.80665;
        
        gFeetWith     = gFeet + gCori;   % walking with rotation (adds outward accel)
        gFeetAgainst  = gFeet - gCori;   % walking against rotation (subtracts)
        
        set(gHUD,'String', sprintf([ ...
            'g_{head}: %.2f g\n' ...
            'g_{feet}: %.2f g\n' ...
            '\\Delta g: %.2f g\n' ...
            'Coriolis @ %.2f m/s (feet): with: %.2f g   against: %.2f g'], ...
            gHead, gFeet, deltaG, WALK_SPEED, gFeetWith, gFeetAgainst));

        
    end

    function onClose(~,~)
        localCleanup(fig, t);
    end

%% ---- slider callbacks ----
    function onRingSlider(src,~)
        val = max(R_OUTER_MIN, min(R_OUTER_MAX, get(src,'Value')));
        set(eRing,'String',sprintf('%.2f',val));
        rebuildStation(val);
    end
    function onRingEdit(src,~)
        val = str2double(get(src,'String'));
        if ~isfinite(val), val = R_outer; end
        val = max(R_OUTER_MIN, min(R_OUTER_MAX, val));
        set(sRing,'Value',val);
        set(src,'String',sprintf('%.2f',val));
        rebuildStation(val);
    end
    function onRpmSlider(src,~)
        val = max(get(src,'Min'), min(get(src,'Max'), get(src,'Value')));
        spinRPM   = val;
        spinRateSec = spinRPM * (2*pi/60);
        set(eRPM,'String',sprintf('%.2f',val));
    end
    function onRpmEdit(src,~)
        val = str2double(get(src,'String'));
        if ~isfinite(val), val = spinRPM; end
        val = max(get(sRPM,'Min'), min(get(sRPM,'Max'), val));
        spinRPM    = val;
        spinRateSec = spinRPM * (2*pi/60);
        set(sRPM,'Value',val);
        set(src,'String',sprintf('%.2f',val));
    end

function onSetRpm1g(~,~)
    % Target RPM to produce 1 g at feet for current R_outer
    g0 = 9.80665;
    R  = max(R_outer, eps);
    targetOmega = sqrt(g0 / R);                  % rad/s
    targetRPM   = targetOmega * 60/(2*pi);       % rev/min

    % Clamp to slider limits
    rpmMin = get(sRPM,'Min'); rpmMax = get(sRPM,'Max');
    rpmNew = min(rpmMax, max(rpmMin, targetRPM));

    spinRPM     = rpmNew;
    spinRateSec = spinRPM * (2*pi/60);
    set(sRPM,'Value',rpmNew);
    set(eRPM,'String',sprintf('%.2f',rpmNew));

    if abs(targetRPM - rpmNew) > 1e-6
        warndlg(sprintf('Required RPM %.2f exceeds slider range [%.2f, %.2f]. Clamped to %.2f.', ...
            targetRPM, rpmMin, rpmMax, rpmNew), 'RPM Clamped', 'nonmodal');
    end
end

function onSetRadius1g(~,~)
    % Target radius to produce 1 g at feet for current omega
    g0 = 9.80665;
    if spinRateSec <= 0
        warndlg('Angular speed is zero; cannot compute radius for 1g. Increase RPM first.', ...
            'Invalid \omega', 'nonmodal');
        return;
    end
    Rneeded = g0 / (spinRateSec^2);              % meters (feet radius)
    Rnew    = min(R_OUTER_MAX, max(R_OUTER_MIN, Rneeded));

    set(sRing,'Value',Rnew);
    set(eRing,'String',sprintf('%.2f',Rnew));
    rebuildStation(Rnew);

    if abs(Rneeded - Rnew) > 1e-6
        warndlg(sprintf(['Required radius %.2f m is outside [%g, %g] m.\n' ...
                         'Clamped to %.2f m.'], ...
            Rneeded, R_OUTER_MIN, R_OUTER_MAX, Rnew), 'Radius Clamped', 'nonmodal');
    end
end

%% ---- geometry rebuild for new ring radius ----
    function rebuildStation(newR)
        % Update globals
        R_outer       = newR;
        R_ring_center = max(R_outer - R_tube, eps);

        % Update torus
        [Xt,Yt,Zt] = torusMesh(R_ring_center, R_tube, U, Vv);
        set(torusSurf,'XData',Xt,'YData',Yt,'ZData',Zt);

        % Recreate spokes (fixed 0.5 m radius)
        if ~isempty(spokeHs), delete(spokeHs(ishandle(spokeHs))); end
        spokeHs = gobjects(0);
        if nSpokes>0
            th = linspace(0,2*pi,nSpokes+1); th(end) = [];
            spokeRad = 0.9;   % <-- fixed radius (meters)
            for ii = 1:nSpokes
                p0 = [0 0 0];
                p1 = [R_ring_center*cos(th(ii)) R_ring_center*sin(th(ii)) 0];
                [SX,SY,SZ] = cylBetween(p0,p1,spokeRad,spokeFacets);
                s = surf(axFG, SX,SY,SZ,'EdgeColor','none', ...
                    'FaceLighting','gouraud', ...
                    'AmbientStrength',0.25,'DiffuseStrength',0.7, ...
                    'SpecularStrength',0.5,'SpecularExponent',25, ...
                    'HitTest','off','PickableParts','none', ...
                    'FaceColor',[0.8 0.8 0.85], 'Parent', G);
                spokeHs(end+1,1) = s; %#ok<AGROW>
            end
        end

        % Recreate hub (depends on centerline R)
        if isgraphics(hubH), delete(hubH); end
        hubHexRadius = 0.25 * R_ring_center;
        [Fhub, Vhub] = regularPrism(6, hubHexRadius, hubDepth);
        hubH = patch('Parent', G, 'Faces', Fhub, 'Vertices', Vhub, ...
            'EdgeColor','none', 'FaceColor', silver, ...
            'FaceLighting','gouraud', ...
            'AmbientStrength',0.25, 'DiffuseStrength',0.7, ...
            'SpecularStrength',0.5, 'SpecularExponent',25, ...
            'HitTest','off','PickableParts','none');

        ZOOM_MIN = R_outer * 0.006;
    end
end

%% ---- helpers ----
function [X,Y,Z] = torusMesh(R, r, U, V)
X = (R + r.*cos(V)) .* cos(U);
Y = (R + r.*cos(V)) .* sin(U);
Z =  r .* sin(V);
end

function [X,Y,Z] = cylBetween(p0,p1,r,n)
if nargin<4, n=24; end
p0 = p0(:)'; p1 = p1(:)';
v = p1 - p0; L = norm(v); if L==0, v=[1 0 0]; L=1; end, v = v / L;
tmp=[0 0 1]; if abs(dot(v,tmp))>0.9, tmp=[0 1 0]; end
u = cross(v,tmp); u = u/norm(u); w = cross(v,u);
theta = linspace(0,2*pi,n+1); theta(end) = [];
circ = r*(cos(theta).' * u + sin(theta).' * w);
P0 = p0 + circ; P1 = p0 + v*L + circ;
X = [P0(:,1) P1(:,1)]; Y = [P0(:,2) P1(:,2)]; Z = [P0(:,3) P1(:,3)];
end

function [F,V] = regularPrism(n, R, depth)
theta = linspace(0, 2*pi, n+1); theta(end) = [];
xy = [R*cos(theta(:)), R*sin(theta(:))];
zTop = +depth/2; zBot = -depth/2;
Vbot = [xy, zBot*ones(n,1)];
Vtop = [xy, zTop*ones(n,1)];
V = [Vbot; Vtop];
Fbot = 1:n;
Ftop = (2*n):-1:(n+1);
Fsides = zeros(n,4);
for i = 1:n
    i2 = mod(i, n) + 1;
    Fsides(i,:) = [i, i2, n+i2, n+i];
end
maxVerts = max(n, 4);
FbotP    = [Fbot, nan(1, maxVerts - numel(Fbot))];
FtopP    = [Ftop, nan(1, maxVerts - numel(Ftop))];
FsidesP  = [Fsides, nan(n, maxVerts - 4)];
F = [FbotP; FtopP; FsidesP];
end

function img = proceduralEarthTexture(h, w)
% Equirectangular "Earth-like" texture with polar ice gradient
    if nargin < 1, h = 1024; end
    if nargin < 2, w = 2048; end

    seed        = 5;
    baseScale   = 24;
    detailScale = 5;
    detailMix   = 0.22;
    sigma_px    = max(h,w)/100;
    landThr     = 0.69;
    thrJitter   = 0.04;

    mergeWin    = 5;
    mergeThr    = 0.58;

    oceanRGB  = uint8([ 20,  80, 160]);
    landRGB   = uint8([ 40, 140,  60]);
    iceRGB    = uint8([230 235 245]);
    iceEdge   = 0.80;
    iceFade   = 0.05;

    rng(seed);
    base   = imresize(rand(ceil(h/baseScale), ceil(w/baseScale)), [h w], 'bilinear');
    detail = imresize(rand(ceil(h/detailScale), ceil(w/detailScale)), [h w], 'bilinear');
    n = (1 - detailMix) * base + detailMix * detail;

    if sigma_px > 0
        x  = -ceil(3*sigma_px):ceil(3*sigma_px);
        g  = exp(-(x.^2)/(2*sigma_px^2)); g = g / sum(g);
        n  = conv2(conv2(n, g, 'same'), g', 'same');
    end

    n = (n - min(n(:))) / (max(n(:)) - min(n(:)) + eps);
    thrMap = landThr + thrJitter * (detail - 0.5);
    land   = n > thrMap;

    if mergeWin > 1
        if mod(mergeWin,2) == 0, mergeWin = mergeWin + 1; end
        kern  = ones(mergeWin) / (mergeWin * mergeWin);
        land  = conv2(double(land), kern, 'same') > mergeThr;
    end

    img = repmat(reshape(oceanRGB,1,1,3), h, w);
    for c = 1:3
        ch  = img(:,:,c);
        ch(land) = [landRGB(c)];
        img(:,:,c) = ch;
    end

    % Polar ice gradient (smoothstep)
    [yy, ~] = ndgrid(linspace(-1,1,h), linspace(-1,1,w));
    t = (abs(yy) - iceEdge) / max(iceFade, eps);
    t = max(0, min(1, t));
    t = t.^2 .* (3 - 2*t);
    t = t .* (0.6 + 0.4*double(land));

    for c = 1:3
        ch = double(img(:,:,c));
        ch = ch.*(1 - t) + double(iceRGB(c)) * t;
        img(:,:,c) = uint8(round(ch));
    end
end

function texOut = bakePlanetLighting(texIn, sunDir, p)
if isempty(texIn), texOut = texIn; return; end
tex = im2double(texIn);  [H,W,~] = size(tex);
sunDir = sunDir(:).'; sunDir = sunDir / norm(sunDir + eps);

lon = linspace(0, 2*pi, W+1); lon(end) = [];
lat = linspace(-pi/2, pi/2, H);
[Lon, Lat] = meshgrid(lon, lat);
Nx = cos(Lat) .* cos(Lon); Ny = cos(Lat) .* sin(Lon); Nz = sin(Lat);
N  = cat(3, Nx, Ny, Nz);
LdotN = max(0, Nx*sunDir(1) + Ny*sunDir(2) + Nz*sunDir(3));
V = cat(3, zeros(H,W), zeros(H,W), ones(H,W));
Hvec(:,:,1) = sunDir(1) + V(:,:,1);
Hvec(:,:,2) = sunDir(2) + V(:,:,2);
Hvec(:,:,3) = sunDir(3) + V(:,:,3);
Hnorm = sqrt(sum(Hvec.^2,3)) + eps;
Hn = bsxfun(@rdivide, Hvec, Hnorm);
HdotN = max(0, sum(Hn .* N, 3));
Spec  = HdotN .^ max(1,p.shiny);
limb = 1 - p.limb*(1 - max(0, Nz).^0.5);
shade = p.ambient + p.diffuse .* LdotN;
shade = shade .* limb;
shadeRGB = repmat(shade, 1, 1, 3) + p.specular * repmat(Spec, 1, 1, 3);
texLit = tex .* shadeRGB;
texLit = max(0, min(1, texLit));
texOut = im2uint8(texLit);
end

function localCleanup(figH, timerH)
if ~isempty(timerH)
    try, if isvalid(timerH), stop(timerH); end, catch, end
    try, if isvalid(timerH), delete(timerH); end, catch, end
end
try
    if isvalid(figH)
        set(figH,'WindowButtonMotionFcn','', ...
                 'WindowButtonDownFcn','', ...
                 'WindowButtonUpFcn','', ...
                 'WindowScrollWheelFcn','');
    end
catch, end
try, if isvalid(figH), delete(figH); end, catch, end
try
    set(0,'ShowHiddenHandles','on');
    delete(findall(0,'Type','figure','Tag','StationSim'));
    oldTimers = timerfindall('Tag','StationSimTimer');
    if ~isempty(oldTimers), stop(oldTimers); delete(oldTimers); end
catch
end
set(0,'ShowHiddenHandles','off');
drawnow;
end

function [gFeet, gHead] = stationGs(R_outer, omega_rad_s, personHeight_m)
g0    = 9.80665;                     
rFeet = max(R_outer, 0);
rHead = max(R_outer - personHeight_m, 0);
gFeet = (omega_rad_s^2 * rFeet) / g0;
gHead = (omega_rad_s^2 * rHead) / g0;
end

function clearOverlay()
    if ~exist('overlayGone','var') || ~overlayGone
        try, if exist('startOverlay','var') && isvalid(startOverlay), delete(startOverlay); end, catch, end
        overlayGone = true; %#ok<NASGU>
    end
end

function tf = pointInRect(pt, rect)
% pt: [x y] in pixels; rect: [x y w h] in pixels
tf = (pt(1) >= rect(1)) && (pt(1) <= rect(1)+rect(3)) && ...
     (pt(2) >= rect(2)) && (pt(2) <= rect(2)+rect(4));
end