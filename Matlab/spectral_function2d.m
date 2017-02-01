function handles = spectral_function2d(handles)

handles.phat = fft2(handles.pd);

    handles.temp = handles.phat .* handles.diffmatrix;



handles.pdiffhat = ifft2(handles.temp);

for i2 = 1 : size(handles.pdiffhat,2)
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
%     handles.udx(i2,:) = handles.udx(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(i2,:)/(3.142*handles.N));
    handles.ud(:,i2) = handles.ud(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(:,i2)/(3.142*handles.N));
end

for i2 = 1 : size(handles.pdiffhat,1)
    if i2 < handles.PMLdepth
        handles.alpha = (1/3)*(((handles.PMLdepth-i2)/ handles.PMLdepth)^3);
    elseif i2 > handles.N - handles.PMLdepth
        handles.alpha = (1/3) * (i2 - ((handles.N-handles.PMLdepth)/handles.PMLdepth)^3);
    else
        handles.alpha = 0;
    end
    handles.ud(i2,:) = handles.ud(i2,:) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(i2,:)/(3.142*handles.N));
%     handles.udy(:,i2) = handles.udy(:,i2) * ((1-handles.alpha)/(1+handles.alpha))-handles.uconst * (1/(1+handles.alpha))*(handles.pdiffhat(:,i2)/(3.142*handles.N));
end

% handles.uhatx = fft2(handles.udx);
% handles.uhaty = fft2(handles.udy);
handles.uhat = fft2(handles.ud);

    handles.temp = handles.uhat .* handles.diffmatrix;


%     temp = uhat .* diffmatrix;
handles.udiffhat = ifft2(handles.temp);
for i2 = 1 : size(handles.udiffhat,1)
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

for i2 = 1 : size(handles.udiffhat,2)
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
%Drive source Locations
handles.pd(ceil(handles.N/2),ceil(handles.N/2)+3) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)+3) +  (1-(handles.src(1,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)+2) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)+2) +  (1-(handles.src(2,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)+1) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)+1) +  (1-(handles.src(3,handles.cntr)));
handles.pd(ceil(handles.N/2)+3,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+3,ceil(handles.N/2)) +  (1-(handles.src(1,handles.cntr)));
handles.pd(ceil(handles.N/2)+2,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+2,ceil(handles.N/2)) +  (1-(handles.src(2,handles.cntr)));
handles.pd(ceil(handles.N/2)+1,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)+1,ceil(handles.N/2)) +  (1-(handles.src(3,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)) +  (1-(handles.src(4,handles.cntr)));
handles.pd(ceil(handles.N/2)-1,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-1,ceil(handles.N/2)) +  (1-(handles.src(5,handles.cntr)));
handles.pd(ceil(handles.N/2)-2,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-2,ceil(handles.N/2)) +  (1-(handles.src(6,handles.cntr)));
handles.pd(ceil(handles.N/2)-3,ceil(handles.N/2)) = handles.pd(ceil(handles.N/2)-3,ceil(handles.N/2)) +  (1-(handles.src(7,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)-1) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)-1) +  (1-(handles.src(5,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)-2) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)-2) +  (1-(handles.src(6,handles.cntr)));
handles.pd(ceil(handles.N/2),ceil(handles.N/2)-3) = handles.pd(ceil(handles.N/2),ceil(handles.N/2)-3) +  (1-(handles.src(7,handles.cntr)));
%Update Counter
handles.cntr = handles.cntr + 1;


end
