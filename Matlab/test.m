clear all;
% close all;
clc;
%prep figs
% figure();
subplot(4,1,1);
%define length
len = 30;
%define diff matrix
tempdiffmatrix = zeros(1,len);
%define input sine
x = sin(-pi:(2*pi)/len:pi);
cosx = cos(-pi:(2*pi)/len:pi);
theor = (diff(x)/((2*pi)/len));
subplot(4,1,1);
plot(x);
hold on;
% plot(theor);
plot(cosx);
legend('x','cosx');
hold off;
for i2 = 1 : len+1
    if i2 <  ceil(len+1/2)
        tempdiffmatrix(i2) =  (i2-1);
    end
    if i2 ==  ceil((len+1)/2)
        tempdiffmatrix(i2) = 0;
    end
    if i2 >  ceil((len+1)/2)
        tempdiffmatrix(i2) = (i2 - (len+1));
    end
end
tempdifft = fft(tempdiffmatrix);
tempdifft = 2*tempdifft/length(tempdifft);
xfft = fft(x);
xfft = 2*xfft/length(xfft);
yfft = real(tempdifft) .* real(xfft);
subplot(4,1,2);
plot(tempdiffmatrix);
hold on;
plot(abs(tempdifft));
plot(abs(xfft));
legend('dempdiffmatrix','tempdifft','xfft');
hold off;
% yfft = xfft .* tempdifft .* (1/len/pi/2);
% for i = 1 : len+1
% yfft(i) = xfft .* tempdiffmatrix;
% end
% yfft = yfft(1:len+1);
subplot(4,1,3);
plot(abs(xfft(1 : round(length(xfft)/2))));
hold on;
plot(abs(yfft(1 : round(length(yfft)/2))));
legend('xfft','yfft');
hold off;
y = ifft(abs(yfft));
y = 2*y/length(y);
y = -(2.*y);
subplot(4,1,4);
plot(x,'--');
hold on;
plot(y);
plot(cosx,':');
legend('x','y','cosx');
hold off;