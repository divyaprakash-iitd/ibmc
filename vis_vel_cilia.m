clear; clc; close all;

% Description: Visualizes velocity fields

% Load u mesh
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');

uFile = dir(strcat('u_0','*'));


pFile = dir(strcat('ib_','*'));
p = load(pFile(1).name);
nc = size(p,1)/4;
nl = 2*nc;

nFiles = length(uFile);

colormap(jet)
v = VideoWriter('ibm.avi','Uncompressed AVI');
open(v);
figure(1)
hold on

for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    contourf(xu,yu,u,50,'edgecolor','none')
    p = load(pFile(iFile).name);
    for i = 1:2:2*nl
        px = p(i,:);
        py = p(i+1,:);
        plot(px,py,'k-o','Markersize',5)
    end
    mesh(xu,yu,0*xu,'FaceAlpha','0.0','EdgeColor','w','LineStyle','-','EdgeAlpha','0.25')
    axis equal
    title(uFile(iFile).name)
    pause(0.001)
    writeVideo(v,getframe(gca));
    if iFile ~= nFiles
        cla
    end
end
close(v)
