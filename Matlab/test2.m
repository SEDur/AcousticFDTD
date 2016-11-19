%S Durbridge, Nov 2016
%Testing spectral differentiation as described in AES conv paper 7963
%GPGPU FDTD By Angus and Caunce
%with information for the differentiator taken from TM/TM/TL Prof. N
%Trefethen
% % http://www.maths.ox.ac.uk/people/nick.trefethen 
% And his wonderful tome
%http://people.maths.ox.ac.uk/trefethen/spectral.html
% 
%This should create a sine wave and cosine wave of a particulat length,
%as well as a frequency domain differentiator of the same length.
%The differentiation matrix is then FFTd to get the impulse response, and
%this is multiplied with the original Sine wave. 
%The output of this is then IFFTd, and should be a cosine that is the same
%as the theoretical cosine

%Init environment
clear all;
close all;
clc;

%prep figs
figure();
subplot(4,3,1);

%%
%First for a sine function

%define length
len = 30;

%define diff matrix
tempdiffmatrix = zeros(1,len);

%define input sine & Analytical result
x = sin((2*pi/len)*(1:len));
cosx = cos((2*pi/len)*(1:len));

%stem the theoretical bits
stem(x);
hold on;
stem(cosx);
legend('x','cosx');
title('Input Signal & Derivative');
hold off;

%Create the differentiator
for i2 = 1 : len
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
%Make the differentiator complex
tempdifft = 1i * tempdiffmatrix;

%FFT the initial sine wave
xfft = fft(x);

%Differentiate the sine wave, by multiplying the sine with the
%differentiator
yfft = tempdifft .* xfft;

%Start stemting
subplot(4,3,4);
%stem the differentiator in the frequ domain
stem(tempdiffmatrix);
title('Noncomplex Differeniator in Frequency Domain');
% hold on;
%stem the time domain idfferentiator
% stem(abs(2 .* tempdifft ./ length(tempdifft)));
%stem the FFT of x
% stem(abs(xfft));
% legend('dempdiffmatrix','tempdifft','xfft');
hold off;
%stem the ffts of X and Y (input and output)
subplot(4,3,7);
stem(abs(xfft(1 : round(length(xfft)/2))));
hold on;
stem(abs(yfft(1 : round(length(yfft)/2))));
legend('xfft','yfft');
title('FFT of original signal and calculated derivative');
hold off;


%Transform y back to time domain
y = ifft(yfft);

%stem final output and compare with theory
subplot(4,3,10);
stem(x,'--');
hold on;
stem(y);
stem(cosx,':');
legend('x','y','cosx');
error = norm(y-cosx,inf);
title(['Result Max Error = ',num2str(error)]);
hold off;

%%
%Now for a hat function
%define length
len = 30;

%define diff matrix
tempdiffmatrix = zeros(1,len);

%define input sine & Analytical result
x = max(0,1-abs(((2*pi/len)*(1:len))-pi)/2);

%stem the theoretical bits
subplot(4,3,2);
stem(x);
hold on;
legend('x');
title('Input Signal & Theoretical Derivative');
hold off;

%Create the differentiator
for i2 = 1 : len
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
%Make the differentiator complex
tempdifft = 1i * tempdiffmatrix;

%FFT the initial sine wave
xfft = fft(x);

%Differentiate the sine wave, by multiplying the sine with the
%differentiator
yfft = tempdifft .* xfft;

%Start stemting
subplot(4,3,5);
%stem the differentiator in the frequ domain
stem(tempdiffmatrix);
title('Noncomplex Differeniator in Frequency Domain');

%stem the ffts of X and Y (input and output)
subplot(4,3,8);
stem(abs(xfft(1 : round(length(xfft)/2))));
hold on;
stem(abs(yfft(1 : round(length(yfft)/2))));
legend('xfft','yfft');
title('FFT of original signal and calculated derivative');
hold off;


%Transform y back to time domain
y = ifft(yfft);

%stem final output and compare with theory
subplot(4,3,11);
stem(x,'--');
hold on;
plot(y);
legend('x','y');
text(2.2,1.4,['max error = ' num2str(error)])
title('Input Vs Output Vs Analytical Res');
hold off;

%%
%Now for a finally for a square wave
%define length
len = 30;

%define diff matrix
tempdiffmatrix = zeros(1,len);

%define input sine & Analytical result
x = zeros(1,len);
x(floor((len/2)-(len/4)):ceil((len/2)+(len/4))) = 1;

%stem the theoretical bits
subplot(4,3,3);
stem(x);
hold on;
legend('x');
title('Input Signal & Theoretical Derivative');
hold off;

%Create the differentiator
for i2 = 1 : len
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
%Make the differentiator complex
tempdifft = 1i * tempdiffmatrix;

%FFT the initial sine wave
xfft = fft(x);

%Differentiate the sine wave, by multiplying the sine with the
%differentiator
yfft = tempdifft .* xfft;

%Start stemting
subplot(4,3,6);
%stem the differentiator in the frequ domain
stem(tempdiffmatrix);
title('Noncomplex Differeniator in Frequency Domain');

%stem the ffts of X and Y (input and output)
subplot(4,3,9);
stem(abs(xfft(1 : round(length(xfft)/2))));
hold on;
stem(abs(yfft(1 : round(length(yfft)/2))));
legend('xfft','yfft');
title('FFT of original signal and calculated derivative');
hold off;


%Transform y back to time domain
y = ifft(yfft);

%stem final output and compare with theory
subplot(4,3,12);
stem(x,'--');
hold on;
plot(y);
legend('x','y');
% text(2.2,1.4,['max error = ' num2str(error)])
title('Input Vs Output Vs Analytical Res');
hold off;

