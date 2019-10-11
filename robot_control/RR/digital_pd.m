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
sampling_time = 0.001;
sampling_int = sampling_time/delta_time;
final_time = 5;
t = 0:delta_time:final_time;
t = t';
t_plot = 0:delta_time:(final_time-2*delta_time);
t_plot = t_plot';
sizes = size(t);

theta1_d = 0.01*sin(2*pi*t/2);
theta2_d = 0.01*cos(2*pi*t/2);
x = L1*cos(theta1_d) + L2*cos(theta1_d+theta2_d);
y = L1*sin(theta1_d) + L2*sin(theta1_d+theta2_d);

%% Calcular a partir dos par�metros desejados os valores de suas derivadas de primeira e segunda ordem

theta1_dotd = (1/delta_time)*diff(theta1_d);
theta2_dotd = (1/delta_time)*diff(theta2_d);
theta1_dot2d = (1/delta_time)*diff(theta1_dotd);
theta2_dot2d = (1/delta_time)*diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

%% Definir condi��es iniciais e ganhos do controlador

q0 = [0.01; 0];
q_dot0 = [0; 0];
gain_p = 100;
gain_v = 2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);
error_plot = zeros(sizes(1)-2,2);
control_plot = zeros(sizes(1)-2,2);
error_old = q_des(:,1) - q0;
error_p_old = q_dot_des(:,1) - q_dot0;
nu=0.1;
%% Laco Principal para Simulacao
for i=1:sizes(1)-2
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    if rem(i-1, sampling_int) == 0
      error = (q_des(:,i) - q0);
      if i==1
        error_p = nu*error_p_old + (error - error_old)*sampling_time;
      else
        error_p = 0;
      end
      control = M*(q_dot2_des(:,i) + k_v.*error_p + k_p.*error) + V + G;
      error_old = error;
      error_p_old = error_p;   
    else
      control = control;
    end
    tau = control;
    error_plot(i,:) = error;
    q_dot2 = inv(M)*(tau-V-G);
    q_dot = q_dot0+ delta_time*q_dot2;
    q = q0 + delta_time*q_dot0 + 0.5*(q_dot2)*(delta_time)^2;
    q0 = q;
    q_dot0 = q_dot;
    control_plot(i,:) = control;
    computed_torque1(i) = control(1);
    computed_torque2(i) = control(2);
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
ylim([-0.01 0.01]);
title('Error- Joints')

