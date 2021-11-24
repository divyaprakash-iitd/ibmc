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


ciliaFile = dir(strcat('ib_loc_c','*'));
particleFile = dir(strcat('ib_loc_p','*'));

cilia = load(ciliaFile(1).name);

ciliaForceFile = dir(strcat('force_ib_loc_c','*'));
particleForceFile = dir(strcat('force_ib_loc_p','*'));

ciliaVelocityFile = dir(strcat('vel_ib_loc_c','*'));
particleVelocityFile = dir(strcat('vel_ib_loc_p','*'));

ncilia = size(cilia,1)/4; % Number of cilia 
% There are two layers in a cilia and thus four rows of values 
% (x1,y1,x2,y2) for all the nodes (equal to the number of columns)
% in the cilia.
nlayers = 2*ncilia; % Total number of layers across all cilia 

% Total number of files (Same for velocity components, cilia and particles)
nFiles = length(uFile);


% Visualize cilia motion over velocity field
colormap(jet)
vid = VideoWriter('ibm.avi','Uncompressed AVI');
open(vid);
figure(1)
fig = gcf;
fig.Position = [1 1 1920 961];
for iFile = 1:nFiles
    %subplot(2,2,1)
    hold on
    
    %% Plot velocity field
    u = load(uFile(iFile).name);
    v = load(vFile(iFile).name);
    % Convert velocity to cell center values
    [uc,vc] = cellcenter(uc,vc,u,v,Nx,Ny);
    umag = sqrt(uc.^2+vc.^2);
    
    contourf(xp,yp,umag,50,'edgecolor','none')
    
    %% Plot cilia
    cilia = load(ciliaFile(iFile).name);
    for i = 1:2:size(cilia,1) % No. of rows in the cilia file
        ciliaX = cilia(i,:); % Read x-coordinates
        ciliaY = cilia(i+1,:); % Read y-coordinates
        plot(ciliaX,ciliaY,'w-o','linewidth',3,'Markersize',5) 
    end
    % Horizontal links
    for i = 1:4:size(cilia,1)
        px1 = cilia(i,:);
        py1 = cilia(i+1,:);
        px2 = cilia(i+2,:);
        py2 = cilia(i+3,:);
        % Horizontal links
        for j = 1:numel(px1)
            plot([px1(j),px2(j)],[py1(j),py2(j)],'w-o','linewidth',1,'Markersize',5) 
        end
        % Diagonal links: 1
        for j = 1:numel(px1)-1
            plot([px1(j),px2(j+1)],[py1(j),py2(j+1)],'w-o','linewidth',1,'Markersize',5)
        end
         % Diagonal links: 2
        for j = 1:numel(px1)-1
            plot([px1(j+1),px2(j)],[py1(j+1),py2(j)],'w-o','linewidth',1,'Markersize',5)
        end
        
    end
    
    %% Plot particles
    particle = load(particleFile(iFile).name);
    for i = 1:2:size(particle,1) % No. of rows in the particle file
        particleX = particle(i,:); % Read x-coordinates
        particleY = particle(i+1,:); % Read y-coordinates
        plot(particleX,particleY,'w-o','linewidth',3,'Markersize',5)
        % Join the ends
        plot([particleX(end), particleX(1)],[particleY(end), particleY(1)],'w-o','linewidth',3,'Markersize',5)
    end
    
   % Horizontal links
    for i = 1:4:size(particle,1) % No. of rows in the particle file
        particleX1 = particle(i,:);
        particleY1 = particle(i+1,:);
        particleX2 = particle(i+2,:);
        particleY2 = particle(i+3,:);
        % Horizontal links
        for j = 1:numel(particleX1)
            plot([particleX1(j),particleX2(j)],[particleY1(j),particleY2(j)],'w-o','linewidth',1,'Markersize',5)
        end
        
        % Diagonal links: 1
        for j = 1:numel(particleX1)-1
            plot([particleX1(j),particleX2(j+1)],[particleY1(j),particleY2(j+1)],'w-o','linewidth',1,'Markersize',5)
        end
        % Join the ends
        plot([particleX1(end),particleX2(1)],[particleY1(end),particleY2(1)],'w-o','linewidth',1,'Markersize',5)
        
         % Diagonal links: 2
        for j = 1:numel(particleX1)-1
            plot([particleX1(j+1),particleX2(j)],[particleY1(j+1),particleY2(j)],'w-o','linewidth',1,'Markersize',5)
        end
        % Join the ends
        plot([particleX1(1),particleX2(end)],[particleY1(1),particleY2(end)],'w-o','linewidth',1,'Markersize',5)
        
    end
    
    %% Plot mesh
    mesh(xp,yp,0*xp,'FaceAlpha','0.0','EdgeColor','w','LineStyle','-','EdgeAlpha','0.25')
    % view(90,0)
    colorbar
    axis equal
    title(uFile(iFile).name)
    
%     %% Plot cilia forces
%     subplot(2,2,2)
%     hold on
%     cilia = load(ciliaFile(iFile).name);
%     ciliaForce = load(ciliaForceFile(iFile).name);
%     for i = 1:2:size(cilia,1) % No. of rows in the cilia file
%         ciliaX = cilia(i,:); % Read x-coordinates
%         ciliaY = cilia(i+1,:); % Read y-coordinates
%         % Read firces
%         ciliaFx = ciliaForce(i,:); % Read x-component force
%         ciliaFy = ciliaForce(i+1,:); % Read y-component force
%         quiver(ciliaX,ciliaY,ciliaFx,ciliaFy,0.25,'k')
%     end
% 
%     %% Plot particle forces
%     particle = load(particleFile(iFile).name);
%     particleForce = load(particleForceFile(iFile).name);
%     for i = 1:2:size(particle,1) % No. of rows in the particle file
%         particleX = particle(i,:); % Read x-coordinates
%         particleY = particle(i+1,:); % Read y-coordinates
%         % Read firces
%         particleFx = particleForce(i,:); % Read x-component force
%         particleFy = particleForce(i+1,:); % Read y-component force
%         quiver(particleX,particleY,particleFx,particleFy,0.25,'k')
%     end
%     
%     %% Plot cilia Velocity
%     subplot(2,2,3)
%     hold on
%     cilia = load(ciliaFile(iFile).name);
%     ciliaVelocity = load(ciliaVelocityFile(iFile).name);
%     for i = 1:2:size(cilia,1) % No. of rows in the cilia file
%         ciliaX = cilia(i,:); % Read x-coordinates
%         ciliaY = cilia(i+1,:); % Read y-coordinates
%         % Read firces
%         ciliaFx = ciliaVelocity(i,:); % Read x-component Velocity
%         ciliaFy = ciliaVelocity(i+1,:); % Read y-component Velocity
%         quiver(ciliaX,ciliaY,ciliaFx,ciliaFy,0.25,'k')
%     end
% 
%     %% Plot particle Velocity
%     particle = load(particleFile(iFile).name);
%     particleVelocity = load(particleVelocityFile(iFile).name);
%     for i = 1:2:size(particle,1) % No. of rows in the particle file
%         particleX = particle(i,:); % Read x-coordinates
%         particleY = particle(i+1,:); % Read y-coordinates
%         % Read firces
%         particleFx = particleVelocity(i,:); % Read x-component Velocity
%         particleFy = particleVelocity(i+1,:); % Read y-component Velocity
%         quiver(particleX,particleY,particleFx,particleFy,0.25,'k')
%     end
    
    title(uFile(iFile).name)
    
%     pause(0.5)
    writeVideo(vid,getframe(gca));
    if iFile ~= nFiles
        clf
    end
end
close(vid)

% Function to convert staggered velocity to cell-centered velocity
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