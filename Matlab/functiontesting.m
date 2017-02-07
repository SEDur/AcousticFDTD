

%%Make The Struct 
%define FS
handles.fs = 1000;

%define density
handles.rho = 1.21;

%define speed of sound
handles.c = 343;

%define total time
handles.T = 10.0;

%define grid width in meters
handles.gridWidth = 10;

%define timestep
handles.dt = 1/(2*handles.fs);
%dfine grid spacing
handles.dx = 2 * handles.dt * handles.c;
%calculate pconst
handles.pconst = handles.rho * handles.c^2 * (handles.dt/handles.dx) * handles.dt * handles.c;
%calculate uconst
handles.uconst = (1/handles.rho)*(handles.dt/handles.dx)*handles.dt*handles.c;
% define pml depth
% handles.PMLdepth = 10;
handles.PMLdepth = 60;
%calc PML depth
% handles.PMLdepth = ceil((abs(handles.gridWidth/handles.dx)/10)/2);
%calc time steps
handles.timestep = abs(handles.T/handles.dt);
%calc grid size
handles.N = ceil(abs(handles.gridWidth/handles.dx)+2*handles.PMLdepth);

%calculate differentiation matrix
handles.tempdiffmatrix = zeros(1,handles.N);
handles.tempdiffmatriy = zeros(handles.N,1);
handles.diffmatrix = zeros(1,handles.N);
handles.diffmatriy = zeros(handles.N,1);
handles.temp = zeros(handles.N, handles.N);
%Calc source
handles.src = ones(7,ceil(handles.T/handles.dt)+10);
% handles.src(1,10:610) = 1 - ((50*10^-10).*sin((2*pi/600)*(1:601)));
% handles.src(2,10:610) = 1 -((60*10^-10).*sin((2*pi/600)*(1:601)));
% handles.src(3,10:610) = 1 -((70*10^-10).*sin((2*pi/600)*(1:601)));
% handles.src(4,10:610) = 1 -((80*10^-10).*sin((2*pi/600)*(1:601)));
% handles.src(5,10:610) = 1 -((70*10^-10).*sin((2*pi/600)*(1:601)));
% handles.src(6,10:610) = 1 -((60*10^-10).*sin((2*pi/600)*(1:601)));
handles.src(7,10:1210) = 1 -((50*10^-10).*sin((2*pi/1200)*(1:1201)));
% handles.src(4,10:110) = ones(101,1);

%Create visualisation data storage matrix
handles.pdstore = zeros(handles.N, handles.N, ceil((handles.T/handles.dt)/100)+1);

%create alpha for PML region
handles.alpha = 0;
%celculate geometry matricies
handles.phat = zeros(handles.N,handles.N);
handles.uhat = zeros(handles.N,handles.N);
handles.pdiffhat = zeros(handles.N,handles.N);
handles.udiffhat = zeros(handles.N,handles.N);
handles.pd = zeros(handles.N,handles.N);
handles.ud = zeros(handles.N,handles.N);
%Create the differentiator
for i2 = 1 : handles.N
    if i2 <  ceil(handles.N+1/2)
        handles.tempdiffmatrix(i2) =  (i2-1);
        handles.tempdiffmatriy(i2,1) =  (i2-1);
    end
    if i2 ==  ceil((handles.N+1)/2)
        handles.tempdiffmatrix(i2) = 0;
        handles.tempdiffmatriy(i2,1) = 0;
    end
    if i2 >  ceil((handles.N+1)/2)
        handles.tempdiffmatrix(i2) = (i2 - (handles.N+1));
        handles.tempdiffmatriy(i2,1) =  (i2 - (handles.N+1));
    end
end

handles.diffmatrix = 1i * handles.tempdiffmatrix;
handles.diffmatriy = 1i * handles.tempdiffmatriy;
handles.hanger = 0;
handles.cntr = 1;

%%
%Run the function



for i = 0 : handles.dt : handles.T
    pause(0.000001) ;
    if handles.hanger == 1
        pause(0.000001) ;
        return
    end
    pause(0.000001) ;
handles = spectral_function(handles);
pause(0.000001) ;
%         if mod(handles.cntr,10)==1
          mesh(real(handles.pd));
%           shading interp;
%         view(2);
        title(sprintf('Time = %.6f s Max = %.6f dB',i ,(20 * log10(max(max(abs(handles.pd))/(10^-12))))));
        drawnow();
%         end
        pause(0.000001) ;
end