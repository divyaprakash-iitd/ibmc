clear; clc; close all;

% Description: Plots the particles of a neighbour's list assigned to each
% cell

cell_data = load('cell.txt');

% Plot the cells
hold on
for i = 1:size(cell_data,1) % No. of rows/entries
    rectangle('Position',cell_data(i,:))
end

% Plot the particles in each cell
p_data = dlmread('cell_particles.txt');
p_data(p_data==0) = nan;

% Plot the particlesplot()
for i = 1:2:size(p_data,1) % No. of rows in the cilia file
    px = p_data(i,:); % Read x-coordinates
    py = p_data(i+1,:); % Read y-coordinates
    plot(px,py,'*','linewidth',3,'Markersize',5) 
end