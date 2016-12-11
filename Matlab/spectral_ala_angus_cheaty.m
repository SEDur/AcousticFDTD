clc;
clear all;
close all;

%define FS
fs = 1000;

%define density
rho = 1.21;

%define speed of sound
c = 343;

%define total time
T = 1.0;

%define grid width in meters
gridWidth = 120;

%define timestep
dt = 1/(2*fs);
%dfine grid spacing
dx = 2 * dt * c;
%calculate pconst
pconst = rho * c^2 * (dt/dx) * dt * c;
%calculate uconst
uconst = (1/rho)*(dt/dx)*dt*c;
% define pml depth 
PMLdepth = 10;
%calc time steps
timestep = abs(T/dt);
%calc grid size
N = ceil(abs(gridWidth/dx)+2*PMLdepth);
%calculate differentiation matrix
tempdiffmatrix = zeros(1,N);
diffmatrix = zeros(N,N);
temp = zeros(1, N);
%Calc source
src = ones(7,ceil(T/dt)+10);
src(1,10:610) = 1 - ((50*10^-10).*sin((2*pi/1200)*(1:601)));
src(2,10:610) = 1 -((70*10^-10).*sin((2*pi/1200)*(1:601)));
src(3,10:610) = 1 -((90*10^-10).*sin((2*pi/1200)*(1:601)));
src(4,10:610) = 1 -((120*10^-10).*sin((2*pi/1200)*(1:601)));
src(5,10:610) = 1 -((90*10^-10).*sin((2*pi/1200)*(1:601)));
src(6,10:610) = 1 -((70*10^-10).*sin((2*pi/1200)*(1:601)));
src(7,10:610) = 1 -((50*10^-10).*sin((2*pi/1200)*(1:601)));

%Create visualisation data storage matrix
pdstore = zeros(ceil(T/dt), N);

%create alpha for PML region
alpha = 0;
%celculate geometry matricies
phat = zeros(1,N);
uhat = zeros(1,N);
pdiffhat = zeros(1,N);
udiffhat = zeros(1,N);
pd = zeros(1,N);
ud = zeros(1,N);
%Create the differentiator
for i2 = 1 : N
    if i2 <  ceil(N+1/2)
        tempdiffmatrix(i2) =  (i2-1);
    end
    if i2 ==  ceil((N+1)/2)
        tempdiffmatrix(i2) = 0;
    end
    if i2 >  ceil((N+1)/2)
        tempdiffmatrix(i2) = (i2 - (N+1));
    end
end

    for i = 1 : length(diffmatrix)
        diffmatrix = 1i * tempdiffmatrix;
    end
figure()
cntr = 1;
cntr2 = 1;
%calculate propagation
for i = 0 : dt : T
    phat = fft(pd);
    temp = phat .* diffmatrix;
    pdiffhat = ifft(temp);
    for i2 = 1 : length(pdiffhat)
        if i2 < PMLdepth
           alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
        elseif i2 > N - PMLdepth
            alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
        else
            alpha = 0;
        end
        ud(i2) = ud(i2) * ((1-alpha)/(1+alpha))-uconst * (1/(1+alpha))*(pdiffhat(i2)/(3.142*N));
        alphastore(cntr) = alpha;
        cntr = cntr + 1;
    end
    
    uhat = fft(ud);
    temp = uhat .* diffmatrix;
    udiffhat = ifft(temp);
    for i2 = 1 : length(udiffhat)
        if i2 < PMLdepth
           alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3); 
        elseif i2 > N - PMLdepth
            alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
        else
            alpha = 0;
        end
        pd(i2) = pd(i2) * ((1-alpha)/(1+alpha))-pconst * (1/(1+alpha))*(udiffhat(i2)/(3.142*N));
        alphastore(cntr) = alpha;
        cntr = cntr + 1;
    end
    
    pd(ceil(N/2)+3) = pd(ceil(N/2)+3) +  (1-(src(1,cntr2)));
    pd(ceil(N/2)+2) = pd(ceil(N/2)+2) +  (1-(src(2,cntr2)));
    pd(ceil(N/2)+1) = pd(ceil(N/2)+1) +  (1-(src(3,cntr2)));
    pd(ceil(N/2)) = pd(ceil(N/2)) +  (1-(src(4,cntr2)));
    pd(ceil(N/2)-1) = pd(ceil(N/2)-1) +  (1-(src(5,cntr2)));
    pd(ceil(N/2)-2) = pd(ceil(N/2)-2) +  (1-(src(6,cntr2)));
    pd(ceil(N/2)-3) = pd(ceil(N/2)-3) +  (1-(src(7,cntr2)));
    subplot(3,1,1:2);
    plot(pd);
    title(sprintf('Time = %.6f s',dt*i));
    subplot(3,1,3);
    plot(alphastore);
    title(sprintf('max alpha = %.3f',max(abs(alphastore))));
    drawnow();
%     pause(0.1);
%     pdstore(cntr,:) = real(pd);
    cntr2 = cntr2 + 1;
end

% for i = 1000 : 100 : length(pdstore)
%     waterfall(pdstore(1:i,100:300)');
%     view([80 30]);
%     title(sprintf('Time = %.6f s',dt*i));
%     drawnow();
% %     pause(0.05);
% end