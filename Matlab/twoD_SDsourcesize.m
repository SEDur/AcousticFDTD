% twoD.m
% Adam Hill
% January 9, 2009

%%Initz Matlab
clear all;

close all;

figure('color', 'w');

%Colormap
CMAP = zeros(256,3);
cl1 = [0 0 1]; %blue
cl2 = [1 1 1]; %white
cl3 = [1 0 0]; %red
 for ncmap = 1 : 128
    fcmap = (ncmap - 1)/128;
    ccmap = (1 - sqrt(fcmap))*cl1 + sqrt(fcmap) * cl2;
    CMAP(ncmap,:)= ccmap;
    ccmap = (1 - fcmap^2)* cl2 + fcmap^2*cl3;
    CMAP(128+ncmap,:) = ccmap;
 end
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
%MAke Movie
makeMovie = 0;
movie_name = 'first2dfdtd.avi';
%Maximum calculation frequency
fmax = 1000 * hertz;
%grid size
gx = (c / fmax) / 10 ;
%Dims
%Dim Size (m)
lx = 10*meters;
ly = 10*meters;

xcells = ceil(lx/gx);
ycells = ceil(lx/gx);

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.4;
alphaR = 0.8;
alphaF = 0.5;
alphaB = 0.25;

%number of sources
snum = 1;
%source locations (x y z))
s1posloc = [ceil(xcells/2.7) ceil(ycells/6.4)];
s1negloc = [(s1posloc(1)+10) (s1posloc(2))];
% s1posloc = [100 101 102 103 104 105 106 107 108 109; 80 80 80 80 80 80 80 80 80 80];
% s1boundloc = [100 101 102 103 104 105 106 107 108 109; 81 81 81 81 81 81 81 81 81 81];
% s1negloc = [100 101 102 103 104 105 106 107 108 109; 82 82 82 82 82 82 82 82 82 82];
%source frequency
s1Freq = 400;

%source phase
s1Phase = 0;

%Source amplitude 
A = 1;

%source size(m) (x y z)
s1size = [1 1 1];

%Time of sim
dt = 1/ (c*sqrt(3/(gx)^2));
% dt = 3.35563e-4;
T = (seconds / fmax) * 40 ;

% generate the source(s) & determine number of time steps needed

tnum = ceil(T/dt);
source1 = zeros(1,tnum);

%generate source grid size

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
            
            t0 = 10;
            t1 = 0 : dt : T;
            phi = -1.6*pi;
            y = A*sin(2*pi*s1Freq*t1 + phi);
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

if makeMovie > 0
    vidObj = VideoWriter(movie_name);
    open(vidObj);
end
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
    p(s1posloc(:,1),s1posloc(:,2)) = p(s1posloc(:,1),s1posloc(:,2)) + source1(n);
    p(s1negloc(:,1),s1negloc(:,2)) = p(s1negloc(:,1),s1negloc(:,2)) + source2(n);
%     p(s1boundloc(:,1),s1boundloc(:,2)) = 0;
%     power(n) = 20*log10(abs(max(p)));
if ~mod(n,2)
    %PLOTTING SECTION
    if makeMovie > 0
    close all;
    end
%         surf(linex, liney, abs(p));
        imagesc(linex, liney, abs(p));
%         imagesc(linex, liney, p);
        shading interp;
        colormap(CMAP);
%         colorbar();
%         caxis([min(source1) max(source1)])
        title(sprintf('Time = %.2f ms',n*dt*1e3),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylabel('Length (meters)', 'Color', [0 0 0]);
        view(2);
%         view([25.6 61.2]);
        drawnow;
        if makeMovie > 0
        %Make movie
        Frame = getframe();
        writeVideo(vidObj, Frame);
        end
    end
end
if makeMovie > 0
close(vidObj);
end