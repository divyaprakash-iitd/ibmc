clear; clc; close all;

% Data for comparison
load Ghia
yG = Ghia.U(:,2);


%% Reynolds Number (100, 400 & 1000)
% To run a case uncomment one of the following three lines

nx = 30; ny = 30; Re = 100; UG = Ghia.U(:,3); 
%nx = 150; ny = 150; Re = 400; UG = Ghia.U(:,4);
%nx = 200; ny = 200; Re = 1000; UG = Ghia.U(:,5);

%% Compare the velocity

% % Load mesh data
% x = load('x.txt');
% y = load('y.txt');
% xm = load('xm.txt');
% ym = load('ym.txt');

imin = 2;
imax = imin+nx-1;
jmin = 2;
jmax = jmin+ny-1;

Lx = 1;
Ly = 1;

% Create Mesh
x(imin:imax+1) = linspace(0,Lx,nx+1);
y(jmin:jmax+1) = linspace(0,Ly,ny+1);
xm(imin:imax)  = 0.5*(x(imin:imax)+x(imin+1:imax+1));
ym(jmin:jmax)  = 0.5*(y(jmin:jmax)+y(jmin+1:jmax+1));

dx = x(3)-x(2);
dy = y(3)-y(2);

x = x(2:end)';
y = y(2:end)';
xm = xm(2:end)';
ym = ym(2:end)';

ux = x;
uy = [ym(1)-dy;ym;ym(end)+dy];

vx = [xm(1)-dx;xm;xm(end)+dx];
vy = y;

[Vx, Vy] = meshgrid(vx,vy);
[Ux, Uy] = meshgrid(ux,uy);
[Px, Py] = meshgrid(xm,ym);


% Load velocity data
uFile = dir(strcat('u_0','*'));
nFiles = length(uFile);

f = figure(1);
xlabel('u')
ylabel('y')
%legend('Location','SouthEast')
title(sprintf('Re = %d, Grid: (%d X %d)',Re,nx,ny))
grid on

for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    xmid = 0.5;
    yq = linspace(yG(1),yG(end),50);
    xq = xmid*ones(size(yq));

    Uq = interp2(Ux,Uy,u,xq,yq);

    % Comparison Plot
    plot(UG,yG,'.','MarkerSize',20,'DisplayName','Ghia et al. (1982)')
    hold on
    plot(Uq,yq,'Parent',f.CurrentAxes)
    title(uFile(iFile).name)
    if iFile ~= nFiles
        hold off
    end
    pause(0.1)
end
hold on
