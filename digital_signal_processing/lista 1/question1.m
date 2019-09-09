clear
clc
close all
N = 12;
n = linspace(0, 2*N-1, 2*N);

figure(1)
clf
M = 4;
stem(n, sin(2*pi*M*n/N));
xlim([0,2*N-1]);
ylim=[-1.1, 1.1];
grid on;
title('Sinal com M=4');


figure(2)
clf;
M = 5;
stem(n, sin(2*pi*M*n/N));
xlim([0,2*N-1]);
ylim=[-1.1, 1.1];
grid on;
title('Sinal com M=5');

figure(3)
clf;
M = 7;
stem(n, sin(2*pi*M*n/N));
xlim([0,2*N-1]);
ylim=[-1.1, 1.1];
grid on;
title('Sinal com M=7');

figure(4)
clf;
M = 10;
stem(n, sin(2*pi*M*n/N));
xlim([0,2*N-1]);
ylim=[-1.1, 1.1];
grid on
title('Sinal com M=10');
