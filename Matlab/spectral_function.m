function handles = spectral_function(handles)
% %define FS
% fs = 1000;
%
% %define density
% rho = 1.21;
%
% %define speed of sound
% c = 343;
%
% %define total time
% T = 10.0;
%
% %define grid width in meters
% gridWidth = 120;
%
% %define timestep
% dt = 1/(2*fs);
% %dfine grid spacing
% dx = 2 * dt * c;
% %calculate pconst
% pconst = rho * c^2 * (dt/dx) * dt * c;
% %calculate uconst
% uconst = (1/rho)*(dt/dx)*dt*c;
% % define pml depth
% PMLdepth = 10;
% %calc time steps
% timestep = abs(T/dt);
% %calc grid size
% N = ceil(abs(gridWidth/dx)+2*PMLdepth);
% %calculate differentiation matrix
% tempdiffmatrix = zeros(1,N);
% diffmatrix = zeros(1,N);
% temp = zeros(N, N);
% %Calc source
% src = ones(7,ceil(T/dt)+10);
% src(1,10:610) = 1 - ((50*10^-10).*sin((2*pi/1200)*(1:601)));
% src(2,10:610) = 1 -((70*10^-10).*sin((2*pi/1200)*(1:601)));
% src(3,10:610) = 1 -((90*10^-10).*sin((2*pi/1200)*(1:601)));
% src(4,10:610) = 1 -((120*10^-10).*sin((2*pi/1200)*(1:601)));
% src(5,10:610) = 1 -((90*10^-10).*sin((2*pi/1200)*(1:601)));
% src(6,10:610) = 1 -((70*10^-10).*sin((2*pi/1200)*(1:601)));
% src(7,10:610) = 1 -((50*10^-10).*sin((2*pi/1200)*(1:601)));
%
% %Create visualisation data storage matrix
% pdstore = zeros(N, N, ceil((T/dt)/100)+1);
%
% %create alpha for PML region
% alpha = 0;
% %celculate geometry matricies
% phat = zeros(N,N);
% uhat = zeros(N,N);
% pdiffhat = zeros(N,N);
% udiffhat = zeros(N,N);
% pd = zeros(N,N);
% ud = zeros(N,N);
% %Create the differentiator
% for i2 = 1 : N
%     if i2 <  ceil(N+1/2)
%         tempdiffmatrix(i2) =  (i2-1);
%     end
%     if i2 ==  ceil((N+1)/2)
%         tempdiffmatrix(i2) = 0;
%     end
%     if i2 >  ceil((N+1)/2)
%         tempdiffmatrix(i2) = (i2 - (N+1));
%     end
% end
%
% for i = 1 : length(diffmatrix)
%     diffmatrix = 1i * tempdiffmatrix;
% end
%
% cntr = 1;

%calculate propagation
% for i = 0 : handles.dt : (handles.T)
handles.phat = fft2(handles.pd);
for i1 = 1 : length(handles.phat)
    handles.tempx(i1,:) = handles.phat(i1,:) .* handles.diffmatrix;
    handles.tempy(:,i1) = handles.phat(:,i1) .* handles.diffmatrix';
end
handles.pdiffhat = ifft2(handles.temp);

for i2 = 1 : length(handles.pdiffhat)
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
    handles.udx(i2,:) = handles.udx(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(i2,:)/(3.142*handles.N));
%     handles.ud(:,i2) = handles.ud(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(:,i2)/(3.142*handles.N));
end

for i2 = 1 + handles.PMLdepth : length(handles.pdiffhat') - handles.PMLdepth
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
%     handles.ud(i2,:) = handles.ud(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(i2,:)/(3.142*handles.N));
    handles.udy(:,i2) = handles.udy(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(:,i2)/(3.142*handles.N));
end

handles.uhatx = fft2(handles.udx);
handles.uhaty = fft2(handles.udy);
for i1 = 1 : length(handles.phat)
    handles.tempx(i1,:) = handles.uhatx(i1,:) .* handles.diffmatrix;
    handles.tempy(:,i1) = handles.uhaty(:,i1) .* handles.diffmatrix';
end
%     temp = uhat .* diffmatrix;
handles.udiffhat = ifft2(handles.temp);
for i2 = 1 : length(handles.udiffhat)
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
    handles.pd(i2,:) = handles.pd(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.pconst * (1/(1+handles.alpha))*(handles.udiffhat(i2,:)/(3.142*handles.N));
%     handles.pd(:,i2) = handles.pd(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.pconst * (1/(1+handles.alpha))*(handles.udiffhat(:,i2)/(3.142*handles.N));
end

for i2 = 1 + handles.PMLdepth : length(handles.udiffhat) - handles.PMLdepth
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
%     handles.pd(i2,:) = handles.pd(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.pconst * (1/(1+handles.alpha))*(handles.udiffhat(i2,:)/(3.142*handles.N));
    handles.pd(:,i2) = handles.pd(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.pconst * (1/(1+handles.alpha))*(handles.udiffhat(:,i2)/(3.142*handles.N));
end

handles.pd(ceil(handles.N/2)+3,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+3,ceil(handles.N/2)) +  (1-(handles.src(1,handles.cntr)));
handles.pd(ceil(handles.N/2)+2,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+2,ceil(handles.N/2)) +  (1-(handles.src(2,handles.cntr)));
handles.pd(ceil(handles.N/2)+1,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+1,ceil(handles.N/2)) +  (1-(handles.src(3,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)) +  (1-(handles.src(4,handles.cntr)));
handles.pd(ceil(handles.N/2)-1,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-1,ceil(handles.N/2)) +  (1-(handles.src(5,handles.cntr)));
handles.pd(ceil(handles.N/2)-2,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-2,ceil(handles.N/2)) +  (1-(handles.src(6,handles.cntr)));
handles.pd(ceil(handles.N/2)-3,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-3,ceil(handles.N/2)) +  (1-(handles.src(7,handles.cntr)));

%         if mod(handles.cntr,100)==1
%           ribbon(real(pd(100:300,100:300)));
%     %     view([80 30]);
%         title(sprintf('Time = %.6f s',dt*i));
%         drawnow();
%         ribbon(real(pd));
%      pdstore(:,:,cntr)  = real(pd);
%         contourf(real(pd));
%         shading interp;
%         surf(real(pd));
%         view([90 55]);
%         title(sprintf('Time = %.6f s',dt*i));
%         drawnow();
%     T/dt*i
%         end
%     pdstore(cntr,:) = real(pd);
handles.cntr = handles.cntr + 1;
% end

% for i = 1000 : 100 : length(pdstore)
%     ribbon(pdstore(100:300,100:300,i));
% %     view([80 30]);
%     title(sprintf('Time = %.6f s',dt*i));
%     drawnow();
%     pause(0.05);
% end

end
