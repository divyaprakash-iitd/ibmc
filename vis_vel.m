clear; clc; close all;

% Description: Visualizes velocity fields

% Load u mesh
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');

uFile = dir(strcat('u_0','*'));
pFile = dir(strcat('ib_','*'));

nFiles = length(uFile);

colormap(jet)

for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    p = load(pFile(iFile).name);
    contourf(xu,yu,u,50,'edgecolor','none')
    hold on
    plot(p(:,1),p(:,2),'w*','Markersize',10)
    title(uFile(iFile).name)
    pause(0.1)
    if iFile ~= nFiles
        cla
    end
end

