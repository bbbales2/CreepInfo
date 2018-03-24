clear all; clc;

data = csvread("/home/bbales2/chuckwalla/CreepInfo/to_plot_matlab.csv", 1, 0);

thickness = data(:, 1);
nsl = data(:, 2);
nsh = data(:, 3);
median = data(:, 4);

fill([thickness; flipud(thickness)], [nsl; flipud(nsh)], 'b');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on;
plot(thickness, median);
hold off;

