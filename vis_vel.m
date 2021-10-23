clear; clc; close all;

% Description: Visualizes velocity fields

% Load u mesh
xu = load('u_x_mesh.txt');
yu = load('u_y_mesh.txt');

uFile = dir(strcat('u_0','*'));
p11File = dir(strcat('1_1_ib_','*'));
p12File = dir(strcat('1_2_ib_','*'));
p21File = dir(strcat('2_1_ib_','*'));
p22File = dir(strcat('2_2_ib_','*'));
p31File = dir(strcat('3_1_ib_','*'));
p32File = dir(strcat('3_2_ib_','*'));

nFiles = length(uFile);

colormap(jet)
v = VideoWriter('ibm.avi','Uncompressed AVI');
open(v);
figure(1)
hold on
for iFile = 1:nFiles
    u = load(uFile(iFile).name);
    p11 = load(p11File(iFile).name);
    p12 = load(p12File(iFile).name);
    p21 = load(p21File(iFile).name);
    p22 = load(p22File(iFile).name);
    p31 = load(p31File(iFile).name);
    p32 = load(p32File(iFile).name);
    contourf(xu,yu,u,50,'edgecolor','none')
%   p = [p ; p(1,:)];
    plot(p11(:,1),p11(:,2),'k-o','Markersize',5)
    plot(p12(:,1),p12(:,2),'k-o','Markersize',5)
    plot(p21(:,1),p21(:,2),'k-o','Markersize',5)
    plot(p22(:,1),p22(:,2),'k-o','Markersize',5)
    plot(p31(:,1),p31(:,2),'k-o','Markersize',5)
    plot(p32(:,1),p32(:,2),'k-o','Markersize',5)
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