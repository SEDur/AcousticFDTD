function handles = spectral_function(handles)

handles.phat = fft(handles.pd,[],1);
for i1 = 1 : size(handles.phat,1)
    handles.temp(i1,:) = handles.phat(i1,:) .* handles.diffmatrix;
end
handles.temp = ifft(handles.temp,[],1);
handles.temp = fft(handles.temp,[],2);
for i1 = 1 : size(handles.phat,2)
%     handles.temp(i1,:) = handles.phat(i1,:) .* handles.diffmatrix;
    handles.temp(:,i1) = handles.phat(:,i1) .* handles.diffmatriy;
end
handles.pdiffhat = ifft(handles.temp,[],2);

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
% handles.uhat = fft2(handles.ud);
% for i1 = 1 : size(handles.uhat,1)
%     handles.temp(i1,:) = handles.uhat(i1,:) .* handles.diffmatrix;
% end
% for i1 = 1 : size(handles.uhat,2)
%     handles.temp(:,i1) = handles.uhat(:,i1) .* handles.diffmatriy;
% end

handles.uhat = fft(handles.ud,[],1);
for i1 = 1 : size(handles.uhat,1)
    handles.temp(i1,:) = handles.uhat(i1,:) .* handles.diffmatrix;
end
handles.temp = ifft(handles.temp,[],1);
handles.temp = fft(handles.temp,[],2);
for i1 = 1 : size(handles.phat,2)
    handles.temp(i1,:) = handles.uhat(i1,:) .* handles.diffmatrix;
%     handles.temp(:,i1) = handles.uhat(:,i1) .* handles.diffmatriy;
end

handles.udiffhat = ifft(handles.temp,[],2);

%     temp = uhat .* diffmatrix;
% handles.udiffhat = ifft2(handles.temp);
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
