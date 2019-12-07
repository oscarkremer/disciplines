%% Limpar console, variaveis e fechar janelas aberta 
clc;
clear; 
close all;
%% Definir massas e comprimentos dos bra�os do Rob� RR

m1 = 1; 
m2 = 1; 
L1 = 1;
L2 = 1;
%% Definir passo e tempo final da simula��o e pegar informa��es de tamanho dos vetores gerados

delta_time = 0.002;
final_time = 10;
t = 0:delta_time:final_time;
t = t';
t_plot = 0:delta_time:(final_time-2*delta_time);
t_plot = t_plot';
sizes = size(t);

%% Gera��o da trajetoria real para calcular os par�metros desejados
%x = 0.3*cos(0.1*t) + 0.2*ones(sizes(1),1);
%y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);

theta1_d = 0.1*sin(2*pi*t/2);
theta2_d = 0.1*cos(2*pi*t/2);
x = L1*cos(theta1_d) + L2*cos(theta1_d+theta2_d);
y = L1*sin(theta1_d) + L2*sin(theta1_d+theta2_d);

theta1_dotd = (1/delta_time)*diff(theta1_d);
theta2_dotd = (1/delta_time)*diff(theta2_d);
theta1_dot2d = (1/delta_time)*diff(theta1_dotd);
theta2_dot2d = (1/delta_time)*diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

q0_gravity = [0.1; 0];
q0 = [0.1; 0];
q_dot0 = [0; 0];
q_dot0_gravity = [0; 0];
gain_p = 100;
gain_v = 2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
g = 9.8;
q_plot = zeros(sizes(1)-2,2);
q_plot_gravity = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);
computed_torque1_gravity = zeros(sizes(1)-2,2);
computed_torque2_gravity = zeros(sizes(1)-2,2);
error_old = (q_des(:,1) - q0);
error_plot = zeros(sizes(1)-2,2);
error_old_gravity = (q_des(:,1) - q0);
error_plot_gravity = zeros(sizes(1)-2,2);



for i=1:sizes(1)-2
    error_new = (q_des(:,i) - q0);
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    M_hat = [1.1*(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) 1.1*m2*L2^2+m2*L1*L2*cos(q0(2)); 
         1.1*m2*L2^2 + m2*L1*L2*cos(q0(2)) 0.8*m2*L2^2]; %erro no termo 1-1 e 2-2 da matriz
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    tau = M_hat*(q_dot2_des(:,i)+k_v.*(q_dot_des(:,i)-q_dot0)+k_p.*(q_des(:,i)-q0))+V+G;
    q_dot2 = inv(M)*(tau-V-G);
    q_dot = q_dot0+ delta_time*q_dot2;
    q = q0 + delta_time*q_dot + 0.5*(q_dot2)*(delta_time)^2;
    error_old = error_new;
    error_plot(i,:) = error_new;
    q0 = q;
    q_dot0 = q_dot;
    computed_torque1(i) = tau(1);
    computed_torque2(i) = tau(2);
    q_plot(i,:) = q;
end

%% Plots para visualiza��o dos �ngulos, posi��o, torque e erros

figure(3)
subplot(2,1,1);
plot(t_plot, q_plot(:,1), '-g')
hold on;
plot(t_plot, q_plot_gravity(:,1),'-b')
plot(t, theta1_d, '-r');
hold off;
title('Trajectory Following - First Rotational Joint');
legend('pd-outerloop', 'pd-gravity', 'desired-theta1');

subplot(2,1,2);
plot(t_plot, q_plot(:,2), '-g')
hold on;
plot(t_plot, q_plot_gravity(:,2), '-b')
plot(t, theta2_d, '-r');
hold off;
title('Trajectory Following - Second Rotational Joint')
legend('pd-outerloop', 'pd-gravity', 'desired-theta2');


figure(4)
subplot(2,1,1);
plot(t_plot, computed_torque1, '-g');
hold on;
plot(t_plot, computed_torque1_gravity, '-r');
hold off;
title('Computed Torque');
legend('pd-outerloop', 'pd-gravity');


subplot(2,1,2);
plot(t_plot, computed_torque2, '-g');
hold on;
plot(t_plot, computed_torque2_gravity, '-r');
hold off;
title('Computed Torque');
legend('pd-outerloop', 'pd-gravity');


figure(5);
subplot(2,1,1);
plot(t_plot, error_plot(:,1), '-g');
hold on;
plot(t_plot, error_plot_gravity(:,1), '-b');
hold off;
title('Error- Joints')
legend('pd-outerloop', 'pd-gravity');



subplot(2,1,2)
plot(t_plot, error_plot(:,2), '-g');
hold on;
plot(t_plot, error_plot_gravity(:,2), '-b');
hold off;
title('Error- Joints')
legend('pd-outerloop', 'pd-gravity');





