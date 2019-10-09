%% Limpar tela, console e variáveis
clc;
clear all;
close all;
%% Declaração sistema e modelo de referência
t = linspace(0,2.5,2500);
u1 = [zeros(250,1); ones(500,1); zeros(1250,1); ones(500,1)];
u2 = [zeros(1000,1); ones(500,1); zeros(500,1); ones(500,1)];
s = tf('s');

G11 = ((2*s+1)/(s^2 + 2*s +1));
G12 = 1/(s^2 + 2*s +1);
G21 = -1/(s^2 + 2*s +1);
G22 = (2*s+1)/(s^2 + 2*s +1);
Wm11 = 10/(s+10);
Wm12 = 0;
Wm21 = 0;
Wm22 = 10/(s+10);


G = [G11 G12; G21 G22];
Wm = [Wm11 Wm12; Wm21 Wm22];
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


%% Declaração matrizes theta

figure(1)
[Y1, Tsim, X1] = lsim(G11,u1,t);
[Y2, Tsim, X2] = lsim(G12,u2,t);
[Ym1, Tsim, Xm1] = lsim(Wm11,u1,t);
plot(Tsim,Y1+Y2,'-r');
hold on;
plot(Tsim,u1,'-b');
plot(Tsim,Ym1,'-g');
legend('Resposta Sistema','Entrada','Referencia')
xlim([0 2.5]);
ylim([-0.5 1.5]);
grid on;
xlabel('Time Seconds');
title('Original System');


figure(2)
[Y1, Tsim, X1] = lsim(G21,u1,t);
[Y2, Tsim, X2] = lsim(G22,u2,t);
[Ym2, Tsim, Xm2] = lsim(Wm22,u2,t);
plot(Tsim,Y1+Y2,'-r');
hold on;
plot(Tsim,u2,'-b');
plot(Tsim,Ym2,'-g');
legend('Resposta Sistema','Entrada','Referencia')
xlim([0 2.5]);
ylim([-0.5 1.5]);
grid on;
xlabel('Time Seconds');
title('Original System');

figure(3)
[Y1, Tsim, X1] = lsim(Y(1,1),u1,t);
[Ym2, Tsim, Xm2] = lsim(Wm11,u1,t);
plot(Tsim,Y1,'-r');
hold on;
plot(Tsim,u1,'-b');
plot(Tsim,Ym2,'-g');
legend('Resposta Sistema','Entrada','Referencia')
xlim([0 2.5]);
ylim([-0.5 1.5]);
grid on;
xlabel('Time Seconds');
title('Original System');

figure(4)
[Y1, Tsim, X1] = lsim(Y(2,2),u2,t);
[Ym2, Tsim, Xm2] = lsim(Wm11,u2,t);
plot(Tsim,Y1,'-r');
hold on;
plot(Tsim,u2,'-b');
plot(Tsim,Ym2,'-g');
legend('Resposta Sistema','Entrada','Referencia')
xlim([0 2.5]);
ylim([-0.5 1.5]);
grid on;
xlabel('Time Seconds');
title('Original System');

