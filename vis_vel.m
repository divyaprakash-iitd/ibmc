clear; clc; close all;

% Description: Visualizes velocity fields
% Description: Visualizes velocity fields

uFile = dir(strcat('u','*'));

nFiles = length(uFile);

colormap(jet)
colorbar

for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    contourf(u,50,'edgecolor','none')
    title(uFile(iFile).name)
    drawnow
end



