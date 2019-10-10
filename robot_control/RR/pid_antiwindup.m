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

delta_time = 0.001;
final_time = 10;
t = 0:delta_time:final_time;
t = t';
t_plot = 0:delta_time:(final_time-2*delta_time);
t_plot = t_plot';
sizes = size(t);

theta1_d = 0.1*sin(2*pi*t/2);
theta2_d = 0.1*cos(2*pi*t/2);
x = L1*cos(theta1_d) + L2*cos(theta1_d+theta2_d);
y = L1*sin(theta1_d) + L2*sin(theta1_d+theta2_d);

%% Calcular a partir dos par�metros desejados os valores de suas derivadas de primeira e segunda ordem

theta1_dotd = (1/delta_time)*diff(theta1_d);
theta2_dotd = (1/delta_time)*diff(theta2_d);
q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';

%% Definir condi��es iniciais e ganhos do controlador

q0 = [0.1; 0];
q_dot0 = [0; 0];
gain_p = 2500;
gain_v = 2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
k_i = [1000; 1000];
g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);
limit_torque = 35;
perturb = 0;
error_old = (q_des(:,1) - q0);
error_plot = zeros(sizes(1)-2,2);
control_plot = zeros(sizes(1)-2,2);
%% Laco Principal para Simulacao
tau_control = zeros(sizes(1)-2,2);
tau = zeros(2,1);
int_error = [0;0];
K_a = 0;
for i=1:sizes(1)-2
    error_new = (q_des(:,i) - q0);
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    control = k_i.*int_error + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0);
    if abs(control(1)) > limit_torque
      int_error(1) = int_error(1) + (q_des(1, i) -q0(1) - K_a*(control(1)-sign(control(1))*limit_torque))*delta_time;
      tau(1) = sign(control(1))*limit_torque;
    else
      int_error(1) = int_error(1) + (q_des(1, i) -q0(1))*delta_time;
      tau(1) = control(1);
    end
    if abs(control(2)) > limit_torque
      int_error(2) = int_error(2) + (q_des(2, i) -q0(2) - K_a*(control(2)-sign(control(2))*limit_torque))*delta_time;
      tau(2) = sign(control(2))*limit_torque;
    else
      int_error(2) = int_error(2) + (q_des(2, i) -q0(2))*delta_time;
      tau(2) = control(2);
    end
    q_dot2 = inv(M)*(tau- V- G -[perturb;perturb]);
    q_dot = q_dot0+ delta_time*q_dot2;
    q = q0 + delta_time*q_dot + 0.5*(q_dot2)*(delta_time)^2;
    error_plot(i,:) = error_new;
    q0 = q;
    q_dot0 = q_dot;
    control_plot(i,:) = control;
    computed_torque1(i) = tau(1);
    computed_torque2(i) = tau(2);
    q_plot(i,:) = q;

end

%% Plots para visualiza��o dos �ngulos, posi��o, torque e erros

figure(3)
subplot(2,1,1);
plot(t_plot, q_plot(:,1), t, theta1_d);
title('Trajectory Following - First Rotational Joint')

subplot(2,1,2);
plot(t_plot, q_plot(:,2), t, theta2_d);
title('Trajectory Following - Second Rotational Joint')


figure(4)

subplot(2,1,1);
plot(t_plot, computed_torque1);
hold on;
plot(t_plot, control_plot(:,1));
hold off;
ylim([0 50]);
title('Computed Torque')

subplot(2,1,2);
plot(t_plot, computed_torque2)
hold on;
plot(t_plot, control_plot(:,2));
hold off;
ylim([0 50]);
title('Computed Torque')


figure(5)
plot(t_plot, error_plot(:,1), t_plot, error_plot(:,2));
ylim([-0.1 0.5]);
title('Error- Joints')

