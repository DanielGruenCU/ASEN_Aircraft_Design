clc;clear;close all;

%% Functions
function [A, A_thick, xu, yu, xl, yl] = naca4_area(code, c, N)
% NACA4_AREA  Cross-sectional area of a NACA 4-digit airfoil by integration.
%
%   [A, A_thick, xu, yu, xl, yl] = naca4_area('2412', c, N)
%
%   Inputs:
%     code : 4-digit NACA designation as a string or number (e.g. '2412')
%     c    : chord length (any units; area is returned in units of c^2)
%     N    : number of points per surface (default 400)
%
%   Outputs:
%     A       : area from the closed upper/lower coordinate loop (shoelace).
%               Includes the camber-line rotation of the thickness.
%     A_thick : area from integrating the thickness distribution, 2*int(yt dx).
%               Camber doesn't change this, so it is a good sanity check.
%     xu,yu   : upper surface coordinates (dimensional)
%     xl,yl   : lower surface coordinates (dimensional)
%
%   Note: uses the closed trailing-edge thickness coefficient (-0.1036).
%   The analytic result for this form is A_thick = 0.6809*t*c^2.

    if nargin < 3, N = 400; end
    if nargin < 2, c = 1;   end

    code = char(string(code));
    m = str2double(code(1))   / 100;   % max camber (fraction of chord)
    p = str2double(code(2))   / 10;    % location of max camber (fraction of chord)
    t = str2double(code(3:4)) / 100;   % max thickness (fraction of chord)

    % Cosine spacing clusters points at the LE, where sqrt(x) is steep
    beta = linspace(0, pi, N);
    x = (1 - cos(beta)) / 2;           % 0 -> 1 (normalized by chord)

    % Thickness distribution
    yt = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x.^2 + 0.2843*x.^3 - 0.1036*x.^4);

    % Camber line and slope
    yc   = zeros(size(x));
    dycdx = zeros(size(x));
    if m > 0 && p > 0
        fwd = x < p;
        yc(fwd)    = m/p^2 * (2*p*x(fwd) - x(fwd).^2);
        dycdx(fwd) = 2*m/p^2 * (p - x(fwd));
        yc(~fwd)    = m/(1-p)^2 * ((1 - 2*p) + 2*p*x(~fwd) - x(~fwd).^2);
        dycdx(~fwd) = 2*m/(1-p)^2 * (p - x(~fwd));
    end
    theta = atan(dycdx);

    % Upper and lower surfaces (thickness applied perpendicular to camber line)
    xu = x - yt.*sin(theta);   yu = yc + yt.*cos(theta);
    xl = x + yt.*sin(theta);   yl = yc - yt.*cos(theta);

    % Dimensionalize
    xu = xu*c; yu = yu*c; xl = xl*c; yl = yl*c;

    % Area 1: closed loop (upper surface LE->TE, then lower surface TE->LE)
    xLoop = [xu, fliplr(xl(1:end-1))];
    yLoop = [yu, fliplr(yl(1:end-1))];
    A = polyarea(xLoop, yLoop);

    % Area 2: integrate thickness distribution directly
    A_thick = trapz(x*c, 2*yt*c);

    A = A / 1.0096; % Adjustment based on CAD comparison
end

function [alpha, cl, cd, cm] = getPolar(airfoil, Re)
% airfoil: lowercase airfoiltools name, e.g. 'naca2412' or 's1223-il'
% Re: one of the precomputed values (50000, 100000, 200000, 500000, 1000000)

    % Names for XFOIL-computed polars look like xf-naca2412-il-500000
    if contains(airfoil, '-il')
        name = sprintf('xf-%s-%d', airfoil, Re);
    else
        name = sprintf('xf-%s-il-%d', airfoil, Re);
    end

    url  = ['http://airfoiltools.com/polar/csv?polar=' name];
    file = [name '.csv'];

    if isfile(file) == false
        websave(file, url);
    end

    % Find the header row (the line starting with "Alpha") instead of
    % hard-coding the number of header lines
    lines = readlines(file);
    hdr   = find(startsWith(lines, 'Alpha'), 1);

    T = readtable(file, 'NumHeaderLines', hdr-1, 'VariableNamingRule', 'preserve');
    alpha = T.Alpha;  cl = T.Cl;  cd = T.Cd;  cm = T.Cm;
end


%% INPUTS:
span = 40; % Wingspan [ft]
airfoil_num = '4412'; % 4 Digit NACA Airfoil
Vland = 100; % Landing Velocity [knots]
alpha_tof = 12; % Takeoff Angle of Attack [deg]
cvt = 0.04;
cht = 0.5;
tail_len = 30; % Length from wing quarter chord to tail quarter chord [ft]


%% OUTPUT GUESS:
Wing_Area_list = 10:10:800; % Wing area (sref) guess [ft^2]

Lift = zeros(length(Wing_Area_list),1);
Weight = zeros(length(Wing_Area_list),1);
CL_ref = zeros(length(Wing_Area_list),1);


%% CONSTANTS (standard sea level):
rho = 1.225; % Air density at sea level [kg/m^3]
g =  9.80665; % Gravitational acceleration [m/s^2]
e = 0.85; % Span efficiency factor
mu = 0.00001789; % Dynamic viscosity [Pa*s]
    
    
%% Input Adjustments
b = span * 0.3048; % span [m]
LT = tail_len * 0.3048; % span [m]

Vland = Vland * 0.5144447; % Convert knots to m/s
Vstall = Vland / 1.3; % Stall Velocity [m/s]
Vtof = 1.1 * Vstall; % Takeoff Velocity [m/s]
    
for i = 1:length(Wing_Area_list)
    Wing_Area = Wing_Area_list(i);
    Sref = Wing_Area * 0.0929;  % area [m^2]
    AR = b^2 / Sref; % aspect ratio
    chord = Sref / b; % chord [m]
    Re = rho * Vstall * chord / mu; % Reynolds Number

    %% Wing Sizing
    % Airfoil Data
    [A, A_thick, xu, yu, xl, yl] = naca4_area(airfoil_num, chord, 400);
    
    AirfoilTools_ReSet = [50000, 100000, 200000, 500000, 1000000];
    [~, idx] = min(abs(log(AirfoilTools_ReSet) - log(Re)));   % nearest on a log scale
    bestRe = AirfoilTools_ReSet(idx);
    
    [alpha, cl, cd, cm] = getPolar(append('naca',airfoil_num), bestRe);
    
    % Get coefficient of lift
    lin = alpha >= -2 & alpha <= 6; % check against a plot
    p = polyfit(alpha(lin), cl(lin), 1);
    a0 = p(1)*180/pi; % per rad
    alpha0 = -p(2)/p(1);
    a = a0 / (1 + a0/(pi*e*AR));
    
    CL = a * deg2rad(alpha_tof - alpha0);
    
    % Add flap contributions (fowler + LE flap)
    del_CLmax = 0.9 * 1.73 * (0.8 * b * 0.3 * chord) / Sref;
    CL = CL + del_CLmax;
    
    %% Weight Sizing
    % Fuselage Weight
    disp_fuse = 77.8; % fuselage displacement volume [ft^3]
    W_fuse = disp_fuse * 28.31685; % fuselage displacement weight [litres -> kgs]
    
    % Wing Weight
    disp_wing = 0.35 * A * b;
    W_wing = disp_wing * 1000; % [litres -> kgs]
    
    % Tail Weight
    S_vtail = cvt * b * Sref / LT;
    S_htail = cht * chord * Sref / LT;
    b_vtail = 10; % [ft]
    b_htail = 15; % [ft]
    c_vtail = S_vtail / (b_vtail * 0.3048);
    c_htail = S_htail / (b_htail * 0.3048);

    [A_vtail, A_thick, xu, yu, xl, yl] = naca4_area(airfoil_num, c_vtail, 400);
    [A_htail, A_thick, xu, yu, xl, yl] = naca4_area(airfoil_num, c_htail, 400);

    disp_vtail = A_vtail * S_vtail * 0.25;
    disp_htail = A_vtail * b_vtail * 0.25;

    W_vtail = disp_vtail * 1000;
    W_htail = disp_htail * 1000;
    
    % Total Weight
    Wtot = W_fuse + W_wing + W_vtail + W_htail;
    
    %% Can plane takeoff
    q_inf = 0.5 * rho * Vstall^2;
    Lift(i) = CL * q_inf * Sref;
    Weight(i) = Wtot * g;
    CL_ref(i) = CL;
end

% Convert to lbf
Lift = Lift * 0.2248089;
Weight = Weight * 0.2248089;

indices = find(Lift - Weight > 0);
[max, position] = max(Lift-Weight);
Weight(position);

disp("Suggested Area: "+Wing_Area_list(position)+" ft^2")
disp("CL at stall: "+CL_ref(position))

figure(1);
hold on;
plot(Wing_Area_list,Weight,LineWidth=1.5)
plot(Wing_Area_list,Lift,LineWidth=1.5)
xline(Wing_Area_list(indices(1)),LineWidth=1,LineStyle="--")
xline(Wing_Area_list(indices(end)),LineWidth=1,LineStyle="--")
ax = gca;
ax.YAxis.Exponent = 0;
legend("weight","lift")
ylabel("Total Weight [lbf]")
xlabel("Planform Wing Area [ft^2]")
hold off;