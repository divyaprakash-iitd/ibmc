clear; clc; close all;

% Description: Visualizes velocity fields

u = load('u.txt');
v = load('v.txt');

%vel = sqrt(u.^2+v.^2);

%contour(u','showtext','on')
%contourf(u')

%u = u(:,1:end-1);

contourf(u',50,'edgecolor','none')
colormap(jet)
colorbar
%contourf(vel','edgecolor','none')
%colorbar
%colormap(jet)
