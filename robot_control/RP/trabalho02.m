%%% Codigo para seguimento de trajetoria de um robo RP %%%
%%% Oscar Schmitt Kremer %%%
%%% Controle de robos - 2019-1 %%%
%%% Ganhos proporcionais, derivativos e integrais %%%


clear 
clc
close all
m = 1;
t = 0:0.002:24;
t = t';
t_plot = 0:0.002:23.996;
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
%y = sin(5*t) + 1.1*ones(sizes(1),1);
%x = cos(5*t);

y = 1*ones(sizes(1),1);
x = 1*sin(5*t);
r_d = (x.^2 + y.^2).^0.5;
theta_d = atan2(y, x);
r_dotd = 0.002*diff(r_d);
theta_dotd = 0.002*diff(theta_d);
r_dot2d = 0.002*diff(r_dotd);
theta_dot2d = 0.002*diff(theta_dotd);

q_des = [theta_d r_d]';
q_dot_des = [theta_dotd r_dotd]';
q_dot2_des = [theta_dot2d r_dot2d]';

q0 = [3.141592/2; 1.5];
q_dot0 = [0; 0];
gain_p = 5000;
gain_v = 0.2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
k_i = [3000; 3000];

g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_force = zeros(sizes(1)-2,2);
computed_torque = zeros(sizes(1)-2,2);

error_old = (q_des(:,1) - q0);
tau = [0;0];
error_plot = zeros(sizes(1)-2,2);
for i=1:sizes(1)-2
    M = [m*q0(2)*q0(2) 0; 
         0 m]; 
    V = [2*m*q0(2)*q_dot0(2)*q_dot0(1);
          -m*q0(2)*(q_dot0(1)^2)];
    G = [m*g*q0(2)*cos(q0(1)); m*g*sin(q0(1))];   
    error_new = (q_des(:,i) - q0);
    int_error = error_old + (error_new - error_old)*0.002;
    q_dot2 = q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*error_new + k_i.*int_error;
    q_dot = q_dot0+ 0.002*q_dot2;
    q = q0 + 0.002*q_dot + 0.5*(q_dot2)*(0.002)^2;
    q0 = q;
    q_dot0 = q_dot;
    tau = M*(q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0)) + V +G;
    computed_torque(i) = tau(1);
    computed_force(i) = tau(2);
    error_old = error_new;
    error_plot(i,:) = error_new;
    q_plot(i,:) = q;
end

figure(3)
subplot(2,1,1);
plot(t_plot, q_plot(:,1), t, theta_d);
title('Trajectory Following - Rotational Joint')

subplot(2,1,2);
plot(t_plot, q_plot(:,2), t, r_d);
title('Trajectory Following - Prismatic Joint')


figure(4)
subplot(2,1,1);
plot(t_plot, computed_torque);
title('Computed Torque')

subplot(2,1,2);
plot(t_plot, computed_force);
title('Computed Force')

x_real  = q_plot(:,2).*cos(q_plot(:,1));
y_real  = q_plot(:,2).*sin(q_plot(:,1));

figure(5)
subplot(2,1,1);
plot(t_plot, x_real);
title('X - Axis')

subplot(2,1,2);
plot(t_plot, y_real);
title('Y - Axis')

figure(6)
subplot(2,1,1);
plot(t_plot, error_plot(:, 1));
title('Error - Rotational')

subplot(2,1,2);
plot(t_plot, error_plot(:, 2));
title('Error - Prismatic')

for i=1:sizes(1)-2
   figure(7)
   hold off
   plot([0 x_real(i)],[0 y_real(i)])
   axis([-2.2 2.2  0 2.5]);
   grid on;
   title('RP Robot Arm');
   hold on
   xlabel('x (m)');
   plot(x(i),y(i),'-x');
   pause(0);
end