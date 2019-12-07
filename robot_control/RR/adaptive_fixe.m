clear 
clc
close all
m1 = 0.8; 
m1hat = 0.85; 
m2 = 2.3; 
m2hat = 2.2; 
L1 = 1;
L2 = 1;
delta_t = 0.001;
final_t = 15;
t = 0:delta_t:final_t;
t = t';
t_plot = 0:delta_t:final_t-2*delta_t;
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
%theta1_d = 0.1*sin(2*pi*t/2);
%theta2_d = 0.1*cos(2*pi*t/2);

%x = cos(theta1_d) + cos(theta1_d + theta2_d);
%y = sin(theta1_d) + sin(theta1_d + theta2_d);
%y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);
%x = 0.3*cos(5*t) + 0.2*ones(sizes(1),1);

%y = 1*ones(sizes(1),1);
%x = 1*sin(5*t);
%cos_theta2 = (x.^2 + y.^2 - L1^2*ones(sizes(1),1) - L2^2*ones(sizes(1),1))/(2*L1*L2);
%sin_theta2 = (1-cos_theta2.^2).^0.5;
theta2_d = sin(t);
%atan2(sin_theta2, cos_theta2);
theta1_d = sin(t);
%atan2(y, x) - atan2(L2*sin(theta2_d), L1*ones(sizes(1),1) + L2*cos(theta2_d));


theta1_dotd = (1/delta_t)*diff(theta1_d);
theta2_dotd = (1/delta_t)*diff(theta2_d);
theta1_dot2d = (1/delta_t)*diff(theta1_dotd);
theta2_dot2d = (1/delta_t)*diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

q0 = [0.1; 0];
q_dot0 = [0; 0];
gain_p = 4;
gain_v = 4;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];

g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);

error_old = (q_des(:,1) - q0);
error_plot = zeros(sizes(1)-2,2);
for i=1:sizes(1)-2
    error_new = (q_des(:,i) - q0);
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    Mhat = [(m1hat+m2hat)*L1^2+m2hat*L2^2+2*m2hat*L1*L2*cos(q0(2)) m2hat*L2^2+m2hat*L1*L2*cos(q0(2)); 
         m2hat*L2^2 + m2hat*L1*L2*cos(q0(2)) m2hat*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    Vhat = [-m2hat*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2hat*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    Ghat = [(m1hat+m2hat)*g*L1*cos(q0(1))+m2hat*g*L2*cos(q0(1)+q0(2)); m2hat*g*L2*cos(q0(1)+q0(2))];   
    tau = Mhat*(q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0)) + Vhat +Ghat;
    q_dot2 = inv(M)*(tau-V-G);
    q_dot = q_dot0+ delta_t*q_dot2;
    q = q0 + delta_t*q_dot + 0.5*(q_dot2)*(delta_t)^2;
    error_old = error_new;
    error_plot(i,:) = error_new;
    q0 = q;
    q_dot0 = q_dot;
    computed_torque1(i) = tau(1);
    computed_torque2(i) = tau(2);
    q_plot(i,:) = q;
end



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

