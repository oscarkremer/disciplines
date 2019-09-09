clear
clc
close all

n = linspace(0, 19, 20);
h = n;
x = [ones(6,1); zeros(14,1)];
h = [n(1:6)'; zeros(14,1)];
y = conv(h,x);


figure(1)
clf;
xlim([0, 39]);
stem(linspace(0, 38, 39), y);
title('Conv h x');
grid on;

x1 = [zeros(5,1); ones(6,1); zeros(14,1)];
h1 = [n(1:6)'; zeros(19,1)];
y1 = conv(x1, h1);

figure(2)
clf;
n = linspace(-5,43,49)';
stem(n, y1);
xlim([-5, 39]);
title('Conv h x');
grid on;

