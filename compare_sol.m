clear; clc; close all;

% Description: Plots the velocity profiles for the fundamental solution of
% stoke's flow

% Load the data
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');
xv = load('v_x_mesh.txt');
yv = load('v_y_mesh.txt');
xp = load('p_x_mesh.txt');
yp = load('p_y_mesh.txt');

% Domain size
Lx = max(xu,[],'all');
Ly = max(yv,[],'all');

uFile = dir(strcat('u_0','*'));
vFile = dir(strcat('v_0','*'));
pFile = dir(strcat('ib_','*'));

nFiles = length(uFile);
u = load(uFile(nFiles).name);
v = load(vFile(nFiles).name);


% Store the values at cell centers
N = 100;
uc = zeros(N,N);
vc = zeros(N,N);

for j = 2:N+1
    for i = 1:N
        vc(i,j-1) = v(i,j); 
    end
end

for j = 1:N
    for i = 2:N+1
        uc(i-1,j) = u(i,j); 
    end
end

umag = sqrt(uc.^2+vc.^2);

figure(1)
contourf(xp,yp,umag,50,'edgecolor','none')
colormap(jet)

% Cobtour plot for analytical solution

mu = 1;
F1 = 0.1;
r = @(x,s) norm(x-s);

U11 = @(x,s) 1/4/pi/mu * (-log(r(x,s)) + ((x(1)-s(1))/r(x,s)).^2);
U12 = @(x,s) 1/4/pi/mu * (((x(1)-s(1))/r(x,s)).*((x(2)-s(2))/r(x,s)));
U21 = @(x,s) 1/4/pi/mu * (((x(1)-s(1))/r(x,s)).*((x(2)-s(2))/r(x,s)));
U22 = @(x,s) 1/4/pi/mu * (-log(r(x,s)) + ((x(2)-s(2))/r(x,s)).^2);

N = 200;
X = linspace(0,Lx,N);
Y = linspace(0,Ly,N);

[x,y] = meshgrid(X,Y);
ux = zeros(size(x));
uy = zeros(size(y));

% Point force location
s = [Lx/2, Ly/2];

for j = 1:N
    for i = 1:N
        xx = [x(i,j), y(i,j)];
        ux(i,j) = F1*U11(xx,s);
        uy(i,j) = F1*U21(xx,s);
    end
end

u = sqrt(ux.^2+uy.^2);

figure(2)
contourf(x,y,u,50,'edgecolor','none')
colormap(jet)

% Compare the two results (Analytical and numerical)

% Location of X-planes
np = 5;
X = linspace(Lx/4,3*Lx/4,np);
ymin = Ly/2; ymax = 0.6; N = 100;
Y = linspace(ymin,ymax,N);

figure(3)
for ip = 1:np
    subplot(1,np,ip)
    hold on
    for j = 1:N
        xx = [X(ip), Y(j)];
        ux = F1*U11(xx,s);
        uy = F1*U21(xx,s);
        u = sqrt(ux^2+uy^2);
        plot(Y(j),ux,'rx')
        plot(Y(j),interp2(xp,yp,uc,X(ip),Y(j)),'kx')
        title(sprintf('x = %.2f',X(ip)))
    end
end


% Plot the velocity due to force calculated numerically
% Middle of the domain

figure(4)
for ip = 1:np
    subplot(1,np,ip)
    plot(Y,interp2(xp,yp,uc,X(ip),Y))
    title(sprintf('x = %.2f',X(ip)))
end
