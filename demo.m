

clear
close all
clc

mozart = double(rgb2gray(imread('mozart.jpg')));

addpath(genpath('C:\Users\alexanders\Desktop\shape-from-shading-master'))

% Define mask ( part of image that contains Object)
[M,N] = size(mozart);
mask = zeros(M,N);

[tilt, slant] = Knill_light_direction(mozart);

niter = 200;
z = Tsi_Shah(mozart,tilt, slant, niter);

figure, surface(abs(z))
