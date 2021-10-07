clear; clc; close all;

% Description: Plots the mesh obtained from the Fortran code

% Load u mesh
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');

z = ones(size(xu));

plot3(xu,yu,z,'ko')
view(0,90)
hold on

% Load v mesh
xv = load('v_x_mesh.txt');
yv = load('v_y_mesh.txt');

z = ones(size(xv));

plot3(xv,yv,z,'k*')

% Load p mesh
xp = load('p_x_mesh.txt');
yp = load('p_y_mesh.txt');

z = ones(size(xp));

plot3(xp,yp,z,'k^')