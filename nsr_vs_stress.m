clear all; clc;

data = csvread("nsr_vs_stress.csv", 1, 1);

thickness = data(:, 1);
nsl = data(:, 2); % 5% quantile
nsh = data(:, 3); % 95% quantile
median = data(:, 4);

fill([thickness; flipud(thickness)], [nsl; flipud(nsh)], 'b');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold on;
plot(thickness, median);
hold off;

