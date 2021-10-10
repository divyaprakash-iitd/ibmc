% clear; clc; close all;

% Description: Plots the velocity profiles for the fundamental solution of
% stoke's flow

nu = 1;

U = @(x) 1/4/pi/nu * [log(1/norm(x)) + x(1)*x(1)/norm(x)^2, x(1)*x(2)/norm(x);
                       x(2)*x(1)/norm(x), log(1/norm(x)) + x(2)*x(2)/norm(x)^2];
                   
                   
        
N = 100;
X = linspace(-2,2,N);
Y = linspace(-2,2,N);

[x,y] = meshgrid(X,Y);
ux = zeros(size(x));
uy = zeros(size(y));

for j = 1:N
    for i = 1:N
        temp = U([x(i,j),y(i,j)]);
%         ux(i,j) = 2*temp(1,1)+2*temp(1,2);
%         uy(i,j) = 2*temp(2,1)+2*temp(2,2);
        ux(i,j) = 2*temp(1,1);
        uy(i,j) = 2*temp(1,2);
    end
end

u = sqrt(ux.^2+uy.^2);

% figure(2)
% contourf(x,y,ux,50,'edgecolor','none')
% colormap(jet)

figure(1)
contourf(x,y,u,50,'edgecolor','none')
colormap(jet)

clear;

xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');
xv = load('v_x_mesh.txt');
yv = load('v_y_mesh.txt');
xp = load('p_x_mesh.txt');
yp = load('p_y_mesh.txt');

uFile = dir(strcat('u_0','*'));
vFile = dir(strcat('v_0','*'));
pFile = dir(strcat('ib_','*'));

nFiles = length(uFile);
u = load(uFile(nFiles).name);
v = load(vFile(nFiles).name);

% Extract u and v for p cells
% N = 100;
% u = u(1:N,2:N+1);
% v = v(2:N+1,1:N);
% xu = xu(1:N,2:N+1);
% xv = xv(2:N+1,1:N);

% Generate the points where you want the values of u and v velocity
uq = interp2(xu,yu,u,xp,yp);
vq = interp2(xv,yv,v,xp,yp);

umag = sqrt(uq.^2+vq.^2);

p = load(pFile(nFiles).name);


figure(2)
hold on
contourf(xp,yp,umag,50,'edgecolor','none')
plot(p(:,1),p(:,2),'w-o','Markersize',5)
title(uFile(nFiles).name)

