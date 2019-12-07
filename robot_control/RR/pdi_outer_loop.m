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
final_time = 5;
t = 0:delta_time:final_time;
t = t';
t_plot = 0:delta_time:(final_time-2*delta_time);
t_plot = t_plot';
sizes = size(t);

%% Gera��o da trajetoria real para calcular os par�metros desejados
%x = 0.3*cos(0.1*t) + 0.2*ones(sizes(1),1);
%y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);


y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);
x = 0.3*cos(5*t) + 0.2*ones(sizes(1),1);

cos_theta2 = (x.^2 + y.^2 - L1^2*ones(sizes(1),1) - L2^2*ones(sizes(1),1))/(2*L1*L2);
sin_theta2 = (1-cos_theta2.^2).^0.5;
theta2_d = atan2(sin_theta2, cos_theta2);
theta1_d = atan2(y, x) - atan2(L2*sin(theta2_d), L1*ones(sizes(1),1) + L2*cos(theta2_d));

%% Calcular a partir dos par�metros desejados os valores de suas derivadas de primeira e segunda ordem


theta1_dotd = (1/delta_time)*diff(theta1_d);
theta2_dotd = (1/delta_time)*diff(theta2_d);
theta1_dot2d = (1/delta_time)*diff(theta1_dotd);
theta2_dot2d = (1/delta_time)*diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

%% Definir condi��es iniciais e ganhos do controlador

q0 = [0.1; 0];
q_dot0 = [0; 0];
gain_p = 2500;
gain_v = 2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
k_i = [1000;1000];
g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);

perturb = 5;
error_old = (q_des(:,1) - q0);
error_plot = zeros(sizes(1)-2,2);
int_error = [0;0];
%% La�o Principal para Simula��o

for i=1:sizes(1)-2
    error_new = (q_des(:,i) - q0);
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    tau = M*(q_dot2_des(:,i)+k_v.*(q_dot_des(:,i)-q_dot0)+k_p.*(q_des(:,i)-q0)+k_i.*int_error)+V+G;
    q_dot2 = inv(M)*(tau-V-G);
    int_error = int_error + (q_des(:, i)-q0)*delta_time;
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
plot(t_plot, q_plot(:,1), t, theta1_d);
title('Trajectory Following - First Rotational Joint')

subplot(2,1,2);
plot(t_plot, q_plot(:,2), t, theta2_d);
title('Trajectory Following - Second Rotational Joint')


figure(4)
subplot(2,1,1);
plot(t_plot, computed_torque1);
title('Computed Torque')

subplot(2,1,2);
plot(t_plot, computed_torque2);
title('Computed Force')


x_real1  = L1*cos(q_plot(:,1));
y_real1  = L1*sin(q_plot(:,1));
x_real2  = L1*cos(q_plot(:,1)) + L2*cos(q_plot(:,1)+q_plot(:,2));
y_real2  = L1*sin(q_plot(:,1)) + L2*sin(q_plot(:,1)+q_plot(:,2));


figure(5)
subplot(2,1,1);
plot(t_plot, x_real2);
title('X - Axis')

subplot(2,1,2);
plot(t_plot, y_real2);
title('Y - Axis')


figure(6)
subplot(2,1,1);
plot(t_plot, error_plot(:,1));
title('Error - Joint 1')

subplot(2,1,2);
plot(t_plot, error_plot(:,2));
title('Error - Joint 2')

%% Rotina para anima��o do rob�
csvwrite ('trajectory.csv', q_plot)
trajectory_x = [];
trajectory_y = [];

for i=1:sizes(1)-2
   figure(7)

   hold off
   plot([0 x_real1(i)],[0 y_real1(i)]);
   

   axis([-2 2  -0.5 2.0]);
   grid on;

   title('RR Robot Arm');
   hold on
   xlabel('x (m)');
   plot(x(i),y(i),'-x');

   plot([x_real1(i) x_real2(i)],[y_real1(i) y_real2(i)]);
   trajectory_x = [trajectory_x; x_real2(i)];
   trajectory_y = [trajectory_y; y_real2(i)];
   plot(trajectory_x, trajectory_y);
   pause(0);
end