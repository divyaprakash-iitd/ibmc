clear; clc; close all;

% Description: Visualizes velocity fields

% Load u mesh
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');

uFile = dir(strcat('u_0','*'));
pFile = dir(strcat('ib_','*'));

nFiles = length(uFile);

colormap(jet)
% v = VideoWriter('ibm.avi');
% open(v);
figure(1)
hold on
for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    p = load(pFile(iFile).name);
    contourf(xu,yu,u,50,'edgecolor','none')
%   p = [p ; p(1,:)];
    plot(p(:,1),p(:,2),'w-o','Markersize',5)
    title(uFile(iFile).name)
    pause(0.001)
%     writeVideo(v,getframe(gca));
    if iFile ~= nFiles
        cla
    end
end
% close(v)

