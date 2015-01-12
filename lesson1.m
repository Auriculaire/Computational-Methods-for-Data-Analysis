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
k = (2 * pi / T) * [0:n/2-1, -n/2: -1]; % Why are they order this way?
                                        % Remember the FFT does the
                                        %    butterfly form.
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
noise = 20;
utn = ut + noise * (randn(1,n) + 1i*randn(1,n)); %note both real and imaginary noise

% Also note noise is added to the frequency domain and NOT the time domain.
% This should wlays be the case.

un = ifft(utn);

subplot(2,1,1), plot(t, u, 'k', t, abs(un), 'm')
subplot(2,1,2), plot(ks,abs(fftshift(ut))/ max(abs(fftshift(ut))), 'k', ...
    ks , abs(fftshift(utn))/ max(abs(fftshift(utn))), 'm')
axis([-25 25 0 1])

noise = 30;
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

% Day 3
% Time-averaging over white noise

% The signal observed
subplot(2,1,1), plot(t,u,'r',t,abs(un),'k')
subplot(2,1,2), plot(ks,abs(fftshift(ut)),'r',ks,abs(fftshift(utn)),'k')

% Consider white noise is assumed random with mean zero
% In the average of the signal over several oberservations, the noise
% contribution should go to zero
noise=30;
avg = zeros(1,n);
realizations = 30;
for j = 1:realizations
    utn = ut + noise * (randn(1,n) + 1i * randn(1,n));
    un = ifft(utn);
    avg = avg + utn;
end
avg = abs(fftshift(avg))/realizations;

plot(ks, abs(fftshift(ut)), 'r', ks, avg, 'k')

% Let's make a movie!
noise=30;
avg = zeros(1,n);
realizations = 100;
for j = 1:realizations
    utn = ut + noise * (randn(1,n) + 1i * randn(1,n));
    un = ifft(utn);
    avg = avg + utn;
    plot(ks, abs(fftshift(ut)), 'r', ks, abs(fftshift(avg))/j, 'k')
    pause(0.2)
end

% Question: Couldn't you do this in the time domain as well?
% Answer: It depends.

% Consider the radar scenario. If you're detecting a moving target, the
% signal in the time domain will move through the domain. This signal
% average over time will go to zero and the target will go undetected.

% However the frequency signal does not change.

T = 60;
n = 512;

t2 = linspace(-T/2, T/2, n+1); t=t2(1:n);
k = (2*pi/T)*[0:n/2-1, -n/2:-1];
ks = fftshift(k);

slice=[0:0.5:10];
[T,S] = meshgrid(t, slice);
[K,S] = meshgrid(k, slice);

% Let's look at the time domain signal of simple sinusoidal signal
U = sech(T-10*sin(S)).*exp(1i*0*T);
subplot(2,1,1)
waterfall(T,S,U), colormap([0 0 0]), view(-15, 70)
% Now in the frequency domain
for j=1:length(slice)
    UT(j, :) = abs(fftshift(fft(U(j,:))));
end    
subplot(2,1,2)
waterfall(fftshift(K), S, UT), colormap([0 0 0]), view(-15, 70)

% Let's shift the frequency
U = sech(T-10*sin(S)).*exp(1i*10*T); % note the change
subplot(2,1,1)
waterfall(T,S,U), colormap([0 0 0]), view(-15, 70)
% Now in the frequency domain
for j=1:length(slice)
    UT(j, :) = abs(fftshift(fft(U(j,:))));
end    
subplot(2,1,2)
waterfall(fftshift(K), S, UT), colormap([0 0 0]), view(-15, 70)
% Note the change in the time domain

% These changes come about as a result of the nonzero frequency number,
% which introduces imaginary components to the time-domain signal. The
% imaginary component can be integrated by taking an absolute value.
waterfall(T,S,abs(U)), colormap([0 0 0]), view(-15, 70)

% Now let's introduce noise
noise = 20;
for j=1:length(slice)
    UTN(j, :) = abs(fftshift(fft(U(j,:)) + noise*(randn(1,n) + 1i*randn(1,n))));
end    
waterfall(fftshift(K), S, UTN), colormap([0 0 0]), view(-15, 70)
% what effect does that have on the time-somain signal.
noise = 20;
for j=1:length(slice)
    UTN(j, :) = fft(U(j,:)) + noise*(randn(1,n) + 1i*randn(1,n));
    UN(j,:) = ifft(UTN(j, :));
end
subplot(2,1,2)
waterfall(fftshift(K), S, abs(fftshift(UTN))), view(-15, 70)
subplot(2,1,1)
waterfall(T, S, abs(UN)), view(-15, 70)