%3D implementation
% S Durbridge 
% 2016

%%Initz Matlab
clear all;
% close all;
% figure(1);
% figure( 'color', 'w');
% figure(2);
% figure( 'color', 'w');
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
p0 = 10^-12;
cstab = sqrt(1/3);
%%
%%Hard Code Variables
%Maximum calculation frequency
fmax = 8000 * hertz;
%grid size
gx = c * (1/fmax) / cstab;
gy = c * (1/fmax) / cstab;
gz = c * (1/fmax) / cstab;
%Dims
%Dim Size (m)
lx = 10*meters;
ly = 10*meters;
lz = 10*meters;

xcells = ceil(lx/gx);
ycells = ceil(ly/gy);
zcells = ceil(lz/gz);

%Boundary Absorption Coefs (0 to 1)
alphaL = 0.5;
alphaR = 0.5;
alphaF = 0.5;
alphaB = 0.5;
alphaT = 0.5;
alphaG = 0.5;

%number of sources
snum = 2;
%source locations
% s1loc = [ceil(xcells/3) ceil(ycells/3)];
% s2loc = [ceil(xcells/1.5) ceil(ycells/1.5)];%source frequency
s1loc = [ceil((ly/gy)/2) ceil((lx/gx)/2) ceil((lz/gz)/4)] ;
s2loc = [ceil((ly/gy)/4) ceil((lx/gx)/2) ceil((lz/gz)/4)];%source frequency
s1Freq = 400;
s2Freq = 400;
%source phase
s1Phase = 0;
s2Phase = 0;
%Source amplitude 
A = 1;

%recieves position
% recieverleftloc = [ceil(lx/gx/2.42) ceil(ly/gx/8) ceil((lz/gz)/2)];
% recieverrightloc = [ceil(lx/gx/2.27) ceil(ly/gx/8) ceil((lz/gz)/2)];
recieverleftloc = [floor((ycells/2) - (0.1/gy)) ceil((xcells/2)-2) ceil(zcells/2)];
recieverrightloc = [ceil((ycells/2) + (0.1/gy)) ceil((xcells/2)+2) ceil(zcells/2)];

%Time of sim
% dt = 1/ (c*sqrt(3/(gx)^2));
dt = 1/ (c*sqrt((1/(gy^2))+(1/(gx^2))+(1/(gz^2))));
% dt = 3.35563e-4;
T = 10*seconds ;

% generate the source(s) & determine number of time steps needed

tnum = ceil(T/dt);
source1 = zeros(1,tnum);
source2 = zeros(1,tnum);
% %           t0 = ceil(T/dt) + 1;
%             t0 = 10;
%             t1 = 0 : dt : T;
%             phi = s1Phase*pi;
%             y = A*sin(2*pi*s1Freq*t1 + phi);
%             gain = linspace(0, 1, ceil(length(y)/10));
%             temp = ones(1, length(y));
%             temp(1 : ceil(length(y)/10)) = gain;
%             y = y.*temp;
%             source1(1, t0 : t0 + ceil(T/dt) - 1) = y;
% source1(1,11:1811) = (sin(0:(pi/1800)*2:(2*pi)))*(p0*10^(100/10));
% source2(1,11:1811) = (sin(0:(pi/1800)*2:(2*pi)))*(p0*10^(100/10));

fc = 0.05;     % Cutoff frequency (normalised 0.5=nyquist)
n0 = 30;        % Initial delay (samples)
sigma=sqrt(2*log(2))/(2*pi*(fc/dt));
n=0:tnum;
source1=exp(-dt^2*(n-n0).^2/(2*sigma^2));
for n = 37 : length(source1)
    if(source1(n) < 0)
       source1(n) = 0; 
    end
end
% source1=filter([1 -1], [1 -0.995], source1); % DC Blocker
% source1 = source1 .* (p0*10^(100/10));
% for n = ceil(tnum/10) : 1 : ceil(tnum/10) + 9 
% source1(n) = source1((n-1) * 2);       
% source2(n) = source2((n-1) * 2);
% end
% %             t0 = ceil(T/dt) + 1;
%             t0 = 10;
%             t1 = 0 : dt : T;
%             phi = s2Phase*pi;
%             y = A*sin(2*pi*s2Freq*t1 + phi);
%             gain = linspace(0, 1, ceil(length(y)/10));
%             temp = ones(1, length(y));
%             temp(1 : ceil(length(y)/10)) = gain;
%             y = y.*temp;
%             source2(1, t0 : t0 + ceil(T/dt) - 1) = y;
        
% initialize the velocity and pressure matrices (matrices are set up in a
% y by x fashion to properly display the 2D space (y = rows, x = columns))
p = zeros(ycells - 1, xcells - 1, zcells - 1);
ux = zeros(ycells - 1, xcells, zcells - 1);
uy = zeros(ycells, xcells - 1, zcells - 1);
uz = zeros(ycells - 1, xcells - 1, zcells);

% set up the multiplication constants for the update equations
uCx = dt/(gx*rho);
uCy = dt/(gy*rho);
uCz = dt/(gz*rho);
pCx = c^2*rho*dt/gx;
pCy = c^2*rho*dt/gy;
pCz = c^2*rho*dt/gz;

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
if alphaT == 0
   alphaT = 1e-016; 
end
if alphaG == 0
   alphaG = 1e-016; 
end
% set the characteristic impedances of the walls
ZR = rho*c*(1 + sqrt(1 - alphaR))/(1 - sqrt(1 - alphaR));
ZL = rho*c*(1 + sqrt(1 - alphaL))/(1 - sqrt(1 - alphaL));
ZF = rho*c*(1 + sqrt(1 - alphaF))/(1 - sqrt(1 - alphaF));
ZB = rho*c*(1 + sqrt(1 - alphaB))/(1 - sqrt(1 - alphaB));
ZT = rho*c*(1 + sqrt(1 - alphaT))/(1 - sqrt(1 - alphaT));
ZG = rho*c*(1 + sqrt(1 - alphaG))/(1 - sqrt(1 - alphaG));
% calulcate the coefficients used for the boundary conditions
Rx = rho*gx/dt;
Ry = rho*gy/dt;
Rz = rho*gz/dt;

% plot vectors
linex = linspace(0, lx - gx, xcells-1);
liney = linspace(0, ly - gy, ycells-1);
linez = linspace(0, lz - gz, zcells-1);

%Initialize recording vectors
leftear = zeros(1,tnum);
rightear = zeros(1,tnum);
meanpstore = zeros(1,tnum);

%Set zsliceloc
zslice = (s1loc(3));

% loop to update the velocities and pressures over the time steps, n
n = 1;
% while or((mean(mean(mean(abs(real(p(:,:,:)))))) > (p0 * 10^(60/10))),(n < (1600)))
while or((meanpstore(n) > (max(meanpstore)-60)),(n < (1600)))
    n = n + 1;
    if mod(n,100)
    (100/tnum)*n
%     10*log10(real(max(max(max(abs(p(:,:,:))))))/p0)
    end
    % update the non-boundary condition nodes for velocity
    ux(:, 2:end-1, :) = ux(:, 2:end-1,:) - uCx*(p(:, 2:end,:) - p(:, 1:end-1, :));
    uy(2:end-1, :, :) = uy(2:end-1, :, :) - uCy*(p(2:end, :, :) - p(1:end-1, :, :));
    uz(:, :, 2:end-1) = uz(:, :, 2:end-1) - uCz*(p(:, :, 2:end) - p(:, :, 1:end-1));
    
%             ux(:, 2:end-1, :) = ux(:, 2:end-1, :) - ...
%             uCx*(p(:, 2:end, :) - p(:, 1:end-1, :));
%         uy(2:end-1, :, :) = uy(2:end-1, :, :) - ...
%             uCy*(p(2:end, :, :) - p(1:end-1, :, :));
%         uz(:, :, 2:end-1) = uz(:, :, 2:end-1) - ...
%             uCz*(p(:, :, 2:end) - p(:, :, 1:end-1));


    % update the velocity at the right wall
    ux(:, end, :) = ((Rx - ZR)/(Rx + ZR))*ux(:, end, :) ...
        + (2/(Rx + ZR))*p(:, end, :);

    %update the velocity at the left wall
    ux(:, 1, :) = ((Rx - ZL)/(Rx + ZL))*ux(:, 1, :) - (2/(Rx + ZL))*p(:, 1, :);

    %update the velocity at the top wall
    uy(end, :, :) = ((Ry - ZT)/(Ry + ZT))*uy(end, :, :) ...
        + (2/(Ry + ZT))*p(end, :, :);

    %update the velocity at the bottom wall
    uy(1, :, :) = ((Ry - ZB)/(Ry + ZB))*uy(1, :, :) - (2/(Ry + ZB))*p(1, :, :);
    
    %update the velocity at the ceiling
    uz(:, :, end) = ((Rz - ZT)/(Rz + ZT))*uz(:, :, end) + ...
    (2/(Rz + ZT))*p(:, :, end);

    %update the velocity at the floor
    uz(:, :, 1) = ((Rz - ZG)/(Rz + ZG))*uz(:, :, 1) - ...
    (2/(Rz + ZG))*p(:, :, 1);

    % update the pressure at all nodes
    p = p - pCx*(ux(:, 2:end, :) - ux(:, 1:end-1, :))...
        - pCy*(uy(2:end, :, :) - uy(1:end-1, :, :))...
        - pCz*(uz(:, :, 2:end) - uz(:, :, 1:end-1));
    % set the pressure at the source location
    % NOTE: source vectors for unused drivers will be zeros
    p(s1loc(1),s1loc(2),s1loc(3)) = p(s1loc(1),s1loc(2),s1loc(3)) - source1(n);
%     p(s2loc(1),s2loc(2)) = p(s2loc(1),s2loc(2)) + -source2(n);
%     power(n) = 20*log10(abs(max(p)));
    leftear(n) = 10*log10(abs(real(p(recieverleftloc(1),recieverleftloc(2),recieverleftloc(3))))/p0);
    rightear(n) = 10*log10(abs(real(p(recieverrightloc(1),recieverrightloc(2),recieverrightloc(3))))/p0);
    %PLOTTING SECTION
    leftear = real(10*log10(leftear/p0));
    rightear = real(10*log10(rightear/p0));
    signal(n) = real(10*log10(source1(n)/p0));
    meanpstore(n) = 10*log10(mean(mean(mean(abs(real(p)))))/p0);
        xslice = ceil(xcells/2); 
        yslice = ycells-1; 
%         zslice = (ceil(zcells-1)/2);
%         slice(abs(p),xslice,yslice,zslice)
        zslice = (ceil(zcells-1)/2);
        slice(abs(p),[xcells-1 ceil(xcells/2)],ycells-1,zslice) 
%         shading interp;
        title(sprintf('Time = %.6f s Max P = %.3f dB',n*dt,10*log10(real(max(max(max(abs(p(:,:,:))))))/p0)),...
            'Color',[0 0 0],'FontSize', 14);
        xlabel('Width (meters)', 'Color', [0 0 0]);
        ylabel('Length (meters)', 'Color', [0 0 0]);
%         figure(2);
%         subplot(2,1,1);
%         plot(dt:dt:n*dt, leftear(1:n));
%         hold on;
%         plot(dt:dt:n*dt, rightear(1:n));
%         hold off;
%         legend('left','right')
%         title((sprintf('Current P recieved by listener = %.3f dB',(rightear(n)))),...
%             'Color',[0 0 0],'FontSize', 14);
%         ylim([0 max(signal)]);
%         subplot(2,1,2);
%         plot(dt:dt:n*dt, signal(1:n));
%         title('whats sent out by the source');
%         ylim([-100 max(signal)]);
%         surf(linex, liney, abs(p));
%         shading interp;
%         title(sprintf('Time = %.6f s',n*dt),...
%             'Color',[0 0 0],'FontSize', 14);
%         xlabel('Width (meters)', 'Color', [0 0 0]);
%         ylabel('Length (meters)', 'Color', [0 0 0]);
%         view(2);
%         view([25.6 61.2]);
        drawnow;
        
end
figure(2);
        subplot(2,1,1);
        plot(dt:dt:n*dt, leftear(1:n));
        hold on;
        plot(dt:dt:n*dt, rightear(1:n));
        plot(dt:dt:n*dt, meanpstore(1:n));
        hold off;
        legend('left','right', 'mean over grid')
        title((sprintf('Current P recieved by listener = %.3f dB & The total sim time was %.6f',(rightear(n)),n*dt)),...
            'Color',[0 0 0],'FontSize', 14);
        ylim([0 max(signal)]);
        subplot(2,1,2);
        plot(dt:dt:n*dt, signal(1:n));
        title('whats sent out by the source');
        ylim([-100 max(signal)]);
% leftear = real(10*log10(leftear/p0));
% rightear = real(10*log10(rightear/p0));
% signal = real(10*log10(source1/p0));
% figure;
% ax = gca;
% plot(leftear);
% hold on;
% plot(rightear);
% plot(signal(1:length(leftear)));
% ax.XTickLabel = [0 : (dt)* 10 : length(leftear)* dt];
% hold off;
