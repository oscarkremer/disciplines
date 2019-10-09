%% Limpar tela, console e variáveis
clc;
clear all;
close all;
pkg load control;
%% Declaração sistema e modelo de referência
t = linspace(0,2.5,1000);
u = [zeros(100,1); ones(200,1); zeros(500,1); ones(200,1)];
s = tf('s');
G = (1/(s^2 + 2*s +1))*[(2*s+1) 1; -1 (2*s+1)];
Wm = (10/(s+10))*eye(2);
%% Declaração matrizes theta
theta11 = [-0.5571 -0.5000;0.5000 -0.5571];
theta12 = [-0.2143 -1.5286;1.5286 -0.2143];
theta13 = [0.1571 -0.8429;0.8429 0.1571];
theta21 = [-0.1857 0;0 -0.1857];
theta22 = [-0.3714 0;0 -0.3714];
theta23 = [-0.1857 0;0 -0.1857];
theta31 = [-3.9714 0;0 -3.9714];

%% Declaração de Y a partir de G e thetas
C0 = 5*eye(2);
Lambda = (s+1)^3;
theta1 = s^2*theta11 + s*theta12 + theta13;
theta2 = s^2*theta21 + s*theta22 + theta23;
theta1 = theta1/Lambda;
theta2 = theta2/Lambda;
Y = G*inv(eye(2) - theta31*G - theta2*G - theta1)*C0;

G11 = (2*s+1)/(s^2 + 2*s +1)
figure(0)
lsim(G11,u,t);
%% Plots das saídas do sistem de referencia e do sistema adaptado
figure(1)
subplot(4,1,1);
step((Y(1,1)))
xlim([0 1]);
ylim([0 1]);
title('Resposta Y-11');

subplot(4,1,2);
step((Y(1,2)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Y-12');

subplot(4,1,3); 
step((Y(2,1)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Y-21');

subplot(4,1,4);
step((Y(2,2)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Y-22');

figure(2)
subplot(4,1,1);
step((Wm(1,1)))
xlim([0 1]);
ylim([0 1]);
title('Resposta Wm-11');

subplot(4,1,2);
step((Wm(1,2)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Wm-12');

subplot(4,1,3); 
step((Wm(2,1)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Wm-21');

subplot(4,1,4);
step((Wm(2,2)));
xlim([0 1]);
ylim([0 1]);
title('Resposta Wm-22');