clear; clc; close all;

% Description: Plots the velocity profiles for the fundamental solution of
% stoke's flow

nu = 100;

U = @(x) 1/4/pi/nu * [log(1/norm(x)) + x(1)*x(1)/(norm(x))^2, x(1)*x(2)/(norm(x))^2;
                       x(2)*x(1)/(norm(x))^2, log(1/norm(x)) + x(2)*x(2)/(norm(x))^2];

        
N = 200;
X = linspace(-3,3,N);
Y = linspace(-3,3,N);

[x,y] = meshgrid(X,Y);
ux = zeros(size(x));
uy = zeros(size(y));

for j = 1:N
    for i = 1:N
        temp = U([x(i,j),y(i,j)]);
%         ux(i,j) = 2*temp(1,1)+2*temp(1,2);
%         uy(i,j) = 2*temp(2,1)+2*temp(2,2);
        ux(i,j) = temp(1,1);
        uy(i,j) = temp(2,1);
    end
end

u = sqrt(ux.^2+uy.^2);

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

% Domain size
Lx = max(xu,[],'all');
Ly = max(yv,[],'all');

uFile = dir(strcat('u_0','*'));
vFile = dir(strcat('v_0','*'));
pFile = dir(strcat('ib_','*'));

% nFiles = length(uFile);
nFiles = 20;
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

figure(2)
contourf(xp,yp,umag,50,'edgecolor','none')
colormap(jet)

% Compare the two results (Analytical and numerical)
nu = 100;

U = @(x) 1/4/pi/nu * [log(1/norm(x)) + x(1)*x(1)/(norm(x))^2, x(1)*x(2)/(norm(x))^2;
                       x(2)*x(1)/(norm(x))^2, log(1/norm(x)) + x(2)*x(2)/norm(x)^2];

% Location of X-planes
np =  5;
X = linspace(-1,1,np);
ymin = 0; ymax = 10; N = 100;
Y = linspace(ymin,ymax,N);

figure(3)
for ip = 1:np
    subplot(1,np,ip)
    hold on
    for j = 1:N
        temp = U([X(ip),Y(j)]);
        ux = temp(1,1);
        uy = temp(2,1);
        u = sqrt(ux^2+uy^2);
        plot(Y(j),ux,'rx')
        title(sprintf('x = %.2f',X(ip)))
    end
end


% Plot the velocity due to force calculated numerically
% Middle of the domain
np = 5;
X = linspace(Lx/4,3*Lx/4,np);

ymin = 0.5; ymax = 0.6; N = 100;
Y = linspace(ymin,ymax,N);

figure(4)
for ip = 1:np
    subplot(1,np,ip)
    plot(Y,interp2(xp,yp,uc,X(ip),Y))
    title(sprintf('x = %.2f',X(ip)))
end
