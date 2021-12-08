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