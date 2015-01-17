% Computational Methods for Data Analysis
% University of Washington
% Homework 1

% Your dog fluffy swallowed a marble. The vet suspects that it has now 
% worked its way into the intestines. Using ultrasound, data is obtained 
% concerning the spatial variations in a small area of the intestines where
% the marble is suspected to be. Unfortunately, fluffy keeps moving and the
% internal fluid movement through the intestines generates highly noisy data.

% Do you want your dog to live? In order to save your dog's life you must 
% located and compute the trajectory of the marble.

% Tabula Rasa
clear all; close all; clc;

L = 15; % Spatial Domain 
nK = 64; % number of Fourier Modes
nT = 20; % number of Time steps

% Let's build an equally spaced grid over each axis of the spatial domain
Grid = linspace(-L, L, nK+1); 
Grid = Grid(1:nK);

% Using the grid along each axis, build a mesh over the spacial domain
[X, Y, Z] = meshgrid(Grid, Grid, Grid);

% Similarly, build a grid over each axis in the frequency domain and
% construct a mesh
k = (2 * pi)/(2 * L) * [0 : (nK/2 - 1), (-nK/2) : -1];
ks = fftshift(k);
[Kx, Ky, Kz] = meshgrid(k, k, k); %rows the columns because matlab is great that way
Kx = fftshift(Kx);
Ky = fftshift(Ky);
Kz = fftshift(Kz);


% Load the data for the problem
filename = 'Testdata.mat';
load(filename)

% Find the freuencies of interest by time-averaging in the frequency domain
avgt(:,:,:) = zeros(nK,nK,nK);
for i = 1:nT
     Un = squeeze(reshape(Undata(i,:), nK, nK, nK));
     Utn = fftn(Un); 
     avgt = avgt + Utn; 
end
avgt = avgt / max(abs(avgt(:)));

% Visually inspect whether a set of frequncies have been resolved
isosurface(Kx, Ky, Kz, fftshift(abs(avgt)), 0.85);
axis([-10 10 -10 10 -10 10]), grid on, drawnow ;

% for j = 1:20
%      close all;
%      j
%      isosurface(Kx, Ky, Kz, fftshift(abs(avgt))/max(abs(avgt(:))), 0.05*j);
%      axis([-10 10 -10 10 -10 10]), grid on, drawnow ;
%      pause(1.0)
% end
 
% Extract the frequencies
[val, index] = max(abs(squeeze(reshape(avgt, nK^3, 1, 1))));
[iKy, iKx, iKz] = ind2sub(size(avgt), index); % Because column-major index is a great idea! /sarcasm
kx = k(iKx); ky = k(iKy); kz = k(iKz); 
% I do not know why the x and y indices are reversed when returned from ind2sub
kx
ky
kz
% Build a Gaussian filter
filter = exp(-((Kx - kx).^2)/7) .* exp(-((Ky - ky).^2)/7) .* exp(-((Kz - kz).^2)/7);

% Initialize coordinate vectors
x = zeros([1, nT]);
y = zeros([1, nT]);
z = zeros([1, nT]);

% Denoise in the spatial domain
Un(:,:,:) = reshape(Undata(1,:), nK, nK, nK);
Utn = fftn(Un);
Ut = fftshift(filter) .* Utn;
U = ifftn(Ut);

% Visually  verify the marble is resolved
% for j= 1:20
%     close all;
%     j
%     isosurface(X, Y, Z, abs(U)/ max(abs(U(:))), 0.05*j)
%     axis([ -20, 20, -20, 20, -20, 20]), grid on, drawnow
%     pause(1)
% end

% Watch a movie
for i = 1:nT
       Un(:,:,:) = reshape(Undata(i,:), nK, nK, nK);
       Utn = fftn(Un);
       Ut = fftshift(filter) .* Utn;
       U = ifftn(Ut);
       isosurface(X, Y, Z, abs(U)/ max(abs(U(:))), 0.8)
       axis([ -15, 15, -15, 15, -15, 15]), grid on, drawnow
       pause(1)
       [val, index] = max(abs(U(:)));
       [iY, iX, iZ] = ind2sub(size(U), index);
       % Again, ind2sub returns iX and iY reversed
       % Column-major indexing folks
       x(i) = Grid(iX);
       y(i) = Grid(iY);
       z(i) = Grid(iZ);
end

plot3(x, y, z)
[x(20),y(20), z(20)]
