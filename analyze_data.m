clear; clc; close all;

% Data for comparison
load Ghia
yG = Ghia.U(:,2);


%% Reynolds Number (100, 400 & 1000)
% To run a case uncomment one of the following three lines

nx = 50; ny = 50; Re = 100; UG = Ghia.U(:,3); 
%nx = 150; ny = 150; Re = 400; UG = Ghia.U(:,4);
%nx = 200; ny = 200; Re = 1000; UG = Ghia.U(:,5);

%% Compare the velocity

% Load mesh data
x = load('x.txt');
y = load('y.txt');
xm = load('xm.txt');
ym = load('ym.txt');

dx = x(3)-x(2);
dy = y(3)-y(2);

x = x(2:end);
y = y(2:end);
xm = xm(2:end);
ym = ym(2:end);

ux = x;
uy = [ym(1)-dy;ym;ym(end)+dy];

vx = [xm(1)-dx;xm;xm(end)+dx];
vy = y;

[Vx, Vy] = meshgrid(vx,vy);
[Ux, Uy] = meshgrid(ux,uy);
[Px, Py] = meshgrid(xm,ym);


% Load velocity data
u = load('u.txt');
v = load('v.txt');

u = u';
v = v';

xmid = 0.5;
yq = linspace(yG(1),yG(end),50);
xq = xmid*ones(size(yq));

Uq = interp2(Ux,Uy,u,xq,yq);

% Comparison Plot
figure(1)
xlabel('u')
ylabel('y')
hold on
plot(UG,yG,'.','MarkerSize',20,'DisplayName','Ghia et al. (1982)')
plot(Uq,yq)
legend('Location','SouthEast')
title(sprintf('Re = %d, Grid: (%d X %d)',Re,nx,ny))
grid on

