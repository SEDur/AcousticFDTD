clc;
clear all;
close all;

%define FS
fs = 10000;
%define density
rho = 1.21;
%define speed of sound
c = 343;
%define total time
T = 1.0;
%define grid width
gridWidth = 2;
%define timestep
dt = 1/(2*fs);
%dfine grid spacing
dx = 2 * dt * c;
%calculate pconst
pconst = rho * c^2 * (dt/dx) * dt * c;
%calculate uconst
uconst = (1/rho)*(dt/dx)*dt*c;
%define pml depth 
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
src = zeros(1,ceil(T/dt)*10);
src(4) = 1;
% for i = 1 : length(src);
%    src(i) =  a;
%    if i < ((ceil(T/dt)/2))
%        a = a+1;
%    end
%    if i >= ((ceil(T/dt)/2)
%        a = a - 1;
%    end
% end
%create alpha for PML region
alpha = 0;
%celculate geometry matricies
phat = zeros(1,N);
uhat = zeros(1,N);
pdiffhat = zeros(1,N);
udiffhat = zeros(1,N);
pd = zeros(1,N);
ud = zeros(1,N);
    for i2 = 1 : N-1
        if i2 <  ceil((N-2)/2)
            tempdiffmatrix(i2) =  (i2-1) * (1+0j);
        end
        if i2 ==  ceil((N-1)/2)
            tempdiffmatrix(i2) = 0 * (1+0j);
        end
        if i2 >  ceil((N-1)/2)
            tempdiffmatrix(i2) = (i2 - 1 - N) * (1+0j);
        end
    end
    ptr1 = ceil(N/2);
    for i = 1 : N
        for i2 = 1:N
            diffmatrix(i2, i) = tempdiffmatrix(ptr1);
            ptr1 = ptr1 + 1;
            if ptr1 > N
                ptr1 = 1;
            end
        end
        ptr1 = ptr1 - 1;
        if ptr1 < 1
            ptr1 = N;
        end
    end
    savediffmat = diffmatrix;
    for i = 1 : length(diffmatrix)
        diffmatrix(i, :) = fft(diffmatrix(i, :));
    end
%calculate source term
% a = 0;
% for i = 1 : ceil((T/dt)/10);
%    src(i) =  a;
%    if i < (((ceil(T/dt)/5)/10))
%        a = a + 0.1;
%    end
%    if i >= (((ceil(T/dt)/5)/10))
%        a = a - 0.1;
%    end
% end
src(10:110) = 1;
cntr = 1;
%calculate propagation
for i = 0 : dt : T
    phat = fft(pd);
%     for i2 = 1 : length(phat)
%         temp = real(phat) .* real(diffmatrix(i2,:));
%     end

    temp = real(phat) .* real(diffmatrix(ceil(N/2),:));
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
    end
    uhat = fft(ud);
%     for i2 = 1 : length(phat)
%         temp(i2) = uhat .* diffmatrix(i2,:);
%     end
        temp = real(uhat) .* real(diffmatrix(ceil(N/2),:));

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
    end
    pd(ceil(N/2)) = pd(ceil(N/2)) +  (src(cntr));
    plot(real(pd));
    title(sprintf('Time = %.6f s',i));
    drawnow();
%     pause(1.0);
    cntr = cntr + 1;
end