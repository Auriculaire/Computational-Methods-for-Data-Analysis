% Computational Methods for Data Analysis
% University of Washington
% Dr. Nathan Kutz
% January 2015

% This script is just a follow along to the Matlab commands used during
% the first weeks lectures.

% Tabula Rasa
clear all; close all; clc;

% Day 1:
% Review of Fourier Transforms


L = 20;  % Domain length
n = 128; % Number of positions sampled

x2 = linspace(-L/2, L/2, n+1); x = x2(1:n); % Domain samples
k = (2*pi/L)*[0:n/2-1 -n/2:-1]; % Frequency samples

u = exp(-x.^2); % Function samples
ut = fft(u);    % Fourier transform of function samples

plot(u)
plot(x, u)
plot(ut)
plot(k, ut)
plot(fftshift(ut))
plot(k, abs(fftshift(ut)))
plot(fftshift(k), abs(fftshift(ut)))

% Derivate computed using fourier transform
% Relationship (in latex):
% \hat{f}^{(n)} = (i*k)^{(n)} * \hat{f} 
u = sech(x);
ut = fft(u);
% Analytical derivatives
ud = -sech(x) .* tanh(x);
u2d = sech(x) - 2*sech(x).^3;


uds=ifft( (1i*k) .* ut ); 
u2ds=ifft( (1i*k).^2 .* ut );

plot(x, ud, 'r', x, uds, 'mo')

% Day 2
% FFT continued
T = 30;
n = 512;
t2 = linspace(-T/2, T/2, n+1); t = t2(1:n);
k = (2 * pi / T) * [0:n/2-1, -n/2: -1];
ks = fftshift(k);

u = sech(t);
ut = fft(u);
subplot(2,1,1), plot(t,u)
subplot(2,1,2), plot(ks,ut)
subplot(2,1,2), plot(ks,fftshift(ut))
subplot(2,1,2), plot(ks,abs(fftshift(ut)))
axis([-25 25 0 1])
subplot(2,1,2), plot(ks,abs(fftshift(ut))/ max(abs(fftshift(ut))))

%adding white noise.
noise = 1;
utn = ut + noise * (randn(1,n) + 1i*rand(1,n)); %note both real and imaginary noise

% Also note noise is added to the frequency domain and NOT the time domain.
% This should wlays be the case.

un = ifft(utn);

subplot(2,1,1), plot(t, u, 'k', t, abs(un), 'm')
subplot(2,1,2), plot(ks,abs(fftshift(ut))/ max(abs(fftshift(ut))), 'k', ...
    ks , abs(fftshift(utn))/ max(abs(fftshift(utn))), 'm')
axis([-25 25 0 1])

noise = 20;
rng(725678)
utn = ut + noise * (randn(1,n) + 1i*rand(1,n)); %note both real and imaginary noise
un = ifft(utn);

subplot(2,1,1), plot(t, abs(un), 'm')
subplot(2,1,2), plot(ks , abs(fftshift(utn))/ max(abs(fftshift(utn))), 'm')
axis([-25 25 0 1])
% note the difficulty in distinguishing the signal of interest from the
% ambient noise.

% One way to recover the signal is using a 'filter'
% There are many many means to filter a signal.
% Basic Gaussian filter used

filter = exp(-k.^2);
plot(ks , abs(fftshift(utn))/ max(abs(fftshift(utn))), 'm', ...
    ks, fftshift(filter), 'b')
axis([-25 25 0 1])

utnf = filter .* utn;
plot(ks , abs(fftshift(utn))/ max(abs(fftshift(utn))), 'm', ...
    ks, fftshift(filter), 'b', ...
    ks , abs(fftshift(utnf))/ max(abs(fftshift(utnf))), 'g')
axis([-25 25 0 1])

unf = ifft(utnf);
plot(t, abs(unf), 'g') % a signal is recovered

% but is the signal just an artifact of the filter?
filter2 = exp(-(k-2).^2);
utnf2 = filter2 .* utn;
unf2 = ifft(utnf2);
plot(t, abs(unf2), 'g') % a signal is recovered