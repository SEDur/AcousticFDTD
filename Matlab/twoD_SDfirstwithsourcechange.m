% twoD.m
% Adam Hill
% January 9, 2009

%%Initz Matlab
clear all;

close all;

figure('color', 'w');

%%Initz Variables

%Units

%Distance
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
%Time
seconds     = 1; 
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

% constants
c     = 343 * meters / seconds; %Speed of sound m/s
rho    = 1.21; %Density of air kg/m^3

%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 1 * kilohertz;
%grid size
gx = (c / fmax) / 10 ;
%Dims
%Dim Size (m)
lx = 1*meters;
ly = 1*meters;

xcells = ceil(lx/gx);
ycells = ceil(lx/gx);

%Boundary Absorption Coefs (0 to 1)
alphaL = 1;
alphaR = 1;
alphaF = 1;
alphaB = 1;

%number of sources
snum = 2;
%source locations
s1loc = [ceil(xcells/3) ceil(ycells/3)];
s2loc = [ceil(xcells/1.5) ceil(ycells/1.5)];%source frequency
s1Freq = 400;
s2Freq = 630;
%source phase
s1Phase = 0;
s2Phase = 90;
%Source amplitude 
A = 1;

%Time of sim
dt = 1/ (c*sqrt(3/(gx)^2));
% dt = 3.35563e-4;
T = (seconds / fmax) * 40 ;

% generate the source(s) & determine number of time steps needed

tnum = ceil(T/dt);
source1 = zeros(1,tnum);
        
%             t0 = ceil(T/dt) + 1;
            t0 = 10;
            t1 = 0 : dt : T;
            phi = s1Phase*pi;
            y = A*sin(2*pi*s1Freq*t1 + phi);
            gain = linspace(0, 1, ceil(length(y)/10));
            temp = ones(1, length(y));
            temp(1 : ceil(length(y)/10)) = gain;
            y = y.*temp;
            source1(1, t0 : t0 + ceil(T/dt) - 1) = y;
        
            source2 = zeros(1,tnum);
        
%             t0 = ceil(T/dt) + 1;
            t0 = 10;
            t1 = 0 : dt : T;
            phi = s2Phase*pi;
            y = A*sin(2*pi*s2Freq*t1 + phi);
            gain = linspace(0, 1, ceil(length(y)/10));
            temp = ones(1, length(y));
            temp(1 : ceil(length(y)/10)) = gain;
            y = y.*temp;
            source2(1, t0 : t0 + ceil(T/dt) - 1) = y;
        
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1);
ux = zeros(ycells - 1, xcells);
uy = zeros(ycells, xcells - 1);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gx*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gx;

% set the wall reflection coefficients
% if alphaX = 0, then slightly adjust to avoid infinite characteristic
% impedance.
if alphaR == 0
   alphaR = 1e-016; 
end
if alphaL == 0
   alphaL = 1e-016; 
end
if alphaF == 0
   alphaF = 1e-016; 
end
if alphaB == 0
   alphaB = 1e-016; 
end
% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZT = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));

% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;
Ry = rho*gx/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gx, ycells-1);

% loop to update the velocities and pressures over the time steps, n
for n = 1 : tnum
    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1) = ux(:, 2:end-1) - uCx*(p(:, 2:end) - p(:, 1:end-1));
    uy(2:end-1, :) = uy(2:end-1, :) - uCy*(p(2:end, :) - p(1:end-1, :));

    % update the velocity at the right wall
    ux(:, end) = ((Rx - ZR)/(Rx + ZR))*ux(:, end) ...
        + (2/(Rx + ZR))*p(:, end);

    %update the velocity at the left wall
    ux(:, 1) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1) - (2/(Rx + ZL))*p(:, 1);

    %update the velocity at the top wall
    uy(end, :) = ((Ry - ZT)/(Ry + ZT))*uy(end, :) ...
        + (2/(Ry + ZT))*p(end, :);

    %update the velocity at the bottom wall
    uy(1, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :) - (2/(Ry + ZB))*p(1, :);

    % update the pressure at all nodes
    p = p - pCx*(ux(:, 2:end) - ux(:, 1:end-1))...
        - pCy*(uy(2:end, :) - uy(1:end-1, :));
    % set the pressure at the source location
    % NOTE: source vectors for unused drivers will be zeros
    p(s1loc(1),s1loc(2)) = p(s1loc(1),s1loc(2)) + source1(n);
    p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) + source2(n);
%     power(n) = 20*log10(abs(max(p)));
    %PLOTTING SECTION
%         surf(linex, liney, abs(p));
        surf(linex, liney, abs(p));
        shading interp;
%         colorbar();
%         caxis([min(source1) max(source1)])
        title(sprintf('Time = %.2f ms',n*dt*1e3),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylabel('Length (meters)', 'Color', [0 0 0]);
        view(2);
%         view([25.6 61.2]);
        drawnow;
        
end
