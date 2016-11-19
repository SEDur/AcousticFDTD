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
T = 10.0;

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
diffmatrix = zeros(1,N);
temp = zeros(N, N);
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
pdstore = zeros(N, N, ceil((T/dt)/100)+1);

%create alpha for PML region
alpha = 0;
%celculate geometry matricies
phat = zeros(N,N);
uhat = zeros(N,N);
pdiffhat = zeros(N,N);
udiffhat = zeros(N,N);
pd = zeros(N,N);
ud = zeros(N,N);
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

cntr = 1;

%calculate propagation
for i = 0 : dt : T
    phat = fft2(pd);
    for i1 = 1 : length(phat)
        temp(i1,:) = phat(i1,:) .* diffmatrix;
        temp(:,i1) = phat(:,i1) .* diffmatrix';
    end
    pdiffhat = ifft2(temp);
    
    for i2 = 1 : length(pdiffhat)
        if i2 < PMLdepth
            alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3);
        elseif i2 > N - PMLdepth
            alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
        else
            alpha = 0;
        end
        ud(i2,:) = ud(i2,:) * ((1-alpha)/(1+alpha))-uconst * (1/(1+alpha))*(pdiffhat(i2,:)/(3.142*N));
        ud(:,i2) = ud(:,i2) * ((1-alpha)/(1+alpha))-uconst * (1/(1+alpha))*(pdiffhat(:,i2)/(3.142*N));
    end
    
    uhat = fft2(ud);
    for i1 = 1 : length(phat)
        temp(i1,:) = uhat(i1,:) .* diffmatrix;
        temp(:,i1) = uhat(:,i1) .* diffmatrix';
    end
    %     temp = uhat .* diffmatrix;
    udiffhat = ifft2(temp);
    for i2 = 1 : length(udiffhat)
        if i2 < PMLdepth
            alpha = (1/3)*(((PMLdepth-i2)/ PMLdepth)^3);
        elseif i2 > N - PMLdepth
            alpha = (1/3) * (i2 - ((N-PMLdepth)/PMLdepth)^3);
        else
            alpha = 0;
        end
        pd(i2,:) = pd(i2,:) * ((1-alpha)/(1+alpha))-pconst * (1/(1+alpha))*(udiffhat(i2,:)/(3.142*N));
        pd(:,i2) = pd(:,i2) * ((1-alpha)/(1+alpha))-pconst * (1/(1+alpha))*(udiffhat(:,i2)/(3.142*N));
    end
    
    pd(ceil(N/2)+3,ceil(N/2)) = pd(ceil(N/2)+3,ceil(N/2)) +  (1-(src(1,cntr)));
    pd(ceil(N/2)+2,ceil(N/2)) = pd(ceil(N/2)+2,ceil(N/2)) +  (1-(src(2,cntr)));
    pd(ceil(N/2)+1,ceil(N/2)) = pd(ceil(N/2)+1,ceil(N/2)) +  (1-(src(3,cntr)));
    pd(ceil(N/2),ceil(N/2)) = pd(ceil(N/2),ceil(N/2)) +  (1-(src(4,cntr)));
    pd(ceil(N/2)-1,ceil(N/2)) = pd(ceil(N/2)-1,ceil(N/2)) +  (1-(src(5,cntr)));
    pd(ceil(N/2)-2,ceil(N/2)) = pd(ceil(N/2)-2,ceil(N/2)) +  (1-(src(6,cntr)));
    pd(ceil(N/2)-3,ceil(N/2)) = pd(ceil(N/2)-3,ceil(N/2)) +  (1-(src(7,cntr)));
    
        if mod(cntr,100)==1
          ribbon(real(pd(100:300,100:300)));
    %     view([80 30]);
        title(sprintf('Time = %.6f s',dt*i));
        drawnow();
    %         ribbon(real(pd));
    %      pdstore(:,:,cntr)  = real(pd);
    %         contourf(real(pd));
    %         shading interp;
    %         surf(real(pd));
    %         view([90 55]);
    %         title(sprintf('Time = %.6f s',dt*i));
    %         drawnow();
    %     T/dt*i
        end
    %     pdstore(cntr,:) = real(pd);
    cntr = cntr + 1;
end

% for i = 1000 : 100 : length(pdstore)
%     ribbon(pdstore(100:300,100:300,i));
% %     view([80 30]);
%     title(sprintf('Time = %.6f s',dt*i));
%     drawnow();
%     pause(0.05);
% end