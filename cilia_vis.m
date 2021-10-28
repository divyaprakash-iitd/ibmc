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

nFiles = length(uFile);
u = load(uFile(nFiles).name);
v = load(vFile(nFiles).name);


% Store the values at cell centers
Nx = size(u,2)-1;
Ny = size(u,1)-2;
uc = zeros(Ny,Nx);
vc = zeros(Ny,Nx);


umag = sqrt(uc.^2+vc.^2);


pFile = dir(strcat('ib_loc','*'));
p = load(pFile(1).name);

pfFile = dir(strcat('force_ib_loc','*'));

nc = size(p,1)/4;
nl = 2*nc;

nFiles = length(uFile);


% Visualize cilia motion over velocity field
colormap(jet)
% vid = VideoWriter('ibm.avi','Uncompressed AVI');
% open(vid);
figure(1)
fig = gcf;
fig.Position = [1 1 1920 961];
for iFile = 1:20:nFiles
    subplot(1,2,1)
    hold on
    u = load(uFile(iFile).name);
    v = load(uFile(iFile).name);
    [uc,vc] = cellcenter(uc,vc,u,v,Nx,Ny);
    umag = sqrt(uc.^2+vc.^2);
    contourf(xp,yp,umag,50,'edgecolor','none')
    p = load(pFile(iFile).name);
    for i = 1:2:2*nl
        px = p(i,:);
        py = p(i+1,:);
        plot(px,py,'k-o','Markersize',5)
    end
    mesh(xp,yp,0*xp,'FaceAlpha','0.0','EdgeColor','w','LineStyle','-','EdgeAlpha','0.25')
%     view(90,0)
    colorbar
    axis equal
    title(uFile(iFile).name)
    
    subplot(1,2,2)
    hold on
    pf = load(pfFile(iFile).name);
    for i = 1:2:2*nl    
        % Read locations
        px = p(i,:);
        py = p(i+1,:);
        % Read forces
        pfx = pf(i,:);
        pfy = pf(i+1,:);
        
        quiver(px,py,pfx,pfy,'k')
    end
    %axis equal
    title(uFile(iFile).name)
    
    pause(0.001)
%     writeVideo(vid,getframe(gca));
    if iFile ~= nFiles
        cla
    end
end
% close(vid)

% % Visulaize cilia forces
% colormap(jet)
% % vid = VideoWriter('ibm.avi','Uncompressed AVI');
% % open(vid);
% figure(2)
% hold on
% 
% for iFile = 1:nFiles
%    
%     p = load(pFile(iFile).name);
%     pf = load(pfFile(iFile).name);
%     for i = 1:2:2*nl
%         % Read locations
%         px = p(i,:);
%         py = p(i+1,:);
%         
%         % Read forces
%         pfx = pf(i,:);
%         pfy = pf(i+1,:);
%         
%         quiver(px,py,pfx,pfy)
%     end
%     axis equal
%     title(uFile(iFile).name)
%     pause(0.001)
% %     writeVideo(vid,getframe(gca));
%     if iFile ~= nFiles
%         cla
%     end
% end
% % close(vid)

function [uc,vc] = cellcenter(uc,vc,u,v,Nx,Ny)
    for j = 2:Nx+1
        for i = 1:Ny
            vc(i,j-1) = v(i,j); 
        end
    end

    for j = 1:Nx
        for i = 2:Ny+1
            uc(i-1,j) = u(i,j); 
        end
    end
end