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

ciliaFile = dir(strcat('ib_loc_c','*'));
nFiles = length(ciliaFile);
cilia = load(ciliaFile(1).name);
ciliaForceFile = dir(strcat('force_ib_loc_c','*'));
ciliaVelocityFile = dir(strcat('vel_ib_loc_c','*'));

ncilia = size(cilia,1)/4; % Number of cilia 
% There are two layers in a cilia and thus four rows of values 
% (x1,y1,x2,y2) for all the nodes (equal to the number of columns)
% in the cilia.
nlayers = 2*ncilia; % Total number of layers across all cilia 


% Visualize cilia motion over velocity field
% colormap(jet)
vid = VideoWriter('ibm.avi','Uncompressed AVI');
open(vid);
figure(1)
fig = gcf;
fig.Position = [1 1 1920 961];
for iFile = 1:nFiles
    subplot(1,2,1)
    hold on
    %% Plot cilia
    cilia = load(ciliaFile(iFile).name);
    for i = 1:2:size(cilia,1) % No. of rows in the cilia file
        ciliaX = cilia(i,:); % Read x-coordinates
        ciliaY = cilia(i+1,:); % Read y-coordinates
        plot(ciliaX,ciliaY,'ko','linewidth',3,'Markersize',5) 
    end
    % Horizontal links
    for i = 1:4:size(cilia,1)
        px1 = cilia(i,:);
        py1 = cilia(i+1,:);
        px2 = cilia(i+2,:);
        py2 = cilia(i+3,:);
        % Horizontal links
        for j = 1:numel(px1)
            plot([px1(j),px2(j)],[py1(j),py2(j)],'ko','linewidth',1,'Markersize',5) 
        end
        % Diagonal links: 1
        for j = 1:numel(px1)-1
            plot([px1(j),px2(j+1)],[py1(j),py2(j+1)],'ko','linewidth',1,'Markersize',5)
        end
         % Diagonal links: 2
        for j = 1:numel(px1)-1
            plot([px1(j+1),px2(j)],[py1(j+1),py2(j)],'ko','linewidth',1,'Markersize',5)
        end
        
    end
    
    %% Plot mesh
    hold on
    mesh(xp,yp,0*xp,'FaceAlpha','0.0','EdgeColor','w','LineStyle','-','EdgeAlpha','0.25')
%     view(90,0)
    %colorbar
    axis equal
    
    %% Plot tip deflection
    subplot(1,2,2)
    tipd = 0.5*(cilia(1,end)+cilia(3,end))-0.25;
    plot(iFile,tipd,'rx')
    % Calculate the deflection according to Euler-Bernoulli's beam theory
    F = 0.0010422;
    L = 0.2;
    E = 62.5;
    I = 6.67e-7;
    EBd = F*L^3/3/E/I;
    hold on
%     plot([0 600],[EBd EBd])
    
% 
%     
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

%     pause(0.5)
    writeVideo(vid,getframe(gca));
    if iFile ~= nFiles
%           cla(fig.Children(1))
            cla(fig.Children(2))
%        cla(fig.Children(3))
    end
end
close(vid)
