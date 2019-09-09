clear 
clc
close all
m1 = 1; 
m2 = 1; 
L1 = 0.5;
L2 = 0.5;
t = 0:0.002:24;
t = t';
t_plot = 0:0.002:23.996;
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);
x = 0.3*cos(5*t) + 0.2*ones(sizes(1),1);

%y = 1*ones(sizes(1),1);
%x = 1*sin(5*t);
cos_theta2 = (x.^2 + y.^2 - L1^2*ones(sizes(1),1) - L2^2*ones(sizes(1),1))/(2*L1*L2);
sin_theta2 = (1-cos_theta2.^2).^0.5;
theta2_d = atan2(sin_theta2, cos_theta2);
theta1_d = atan2(y, x) - atan2(L2*sin(theta2_d), L1*ones(sizes(1),1) + L2*cos(theta2_d));
theta1_dotd = diff(theta1_d);
theta2_dotd = diff(theta2_d);
theta1_dot2d = diff(theta1_dotd);
theta2_dot2d = diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

q0 = [3.141592/2; 0];
q_dot0 = [0; 0];
gain_p = 5000;
gain_v = 0.2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
k_i = [3000; 3000];

g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);

error_old = (q_des(:,1) - q0);
tau = [0;0];
error_plot = zeros(sizes(1)-2,2);
for i=1:sizes(1)-2
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    error_new = (q_des(:,i) - q0);
    int_error = error_old + (error_new - error_old)*0.002;
    q_dot2 = q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*error_new + k_i.*int_error;
    q_dot = q_dot0+ 0.002*q_dot2;
    q = q0 + 0.002*q_dot + 0.5*(q_dot2)*(0.002)^2;
    q0 = q;
    q_dot0 = q_dot;
    tau = M*(q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0)) + V +G;
    computed_torque1(i) = tau(1);
    computed_torque2(i) = tau(2);
    error_old = error_new;
    error_plot(i,:) = error_new;
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
plot(t_plot, error_plot(:, 1));
title('Error - Rotational')

subplot(2,1,2);
plot(t_plot, error_plot(:, 2));
title('Error - Prismatic')

figure(7)
subplot(2,1,1);
plot(t_plot, (x_real1.^2 + y_real1.^2).^0.5);
title('Distancia x1-origem')

subplot(2,1,2);
plot(t_plot, ((x_real2-x_real1).^2 + (y_real2-y_real1).^2).^0.5);
title('Distancia x2-x1')


for i=1:sizes(1)-2
   figure(8)
   hold off
   plot([0 x_real1(i)],[0 y_real1(i)]);
   axis([-0.75 0.75  -0.5 1.0]);
   grid on;

   title('RP Robot Arm');
   hold on
   xlabel('x (m)');
   plot([x_real1(i) x_real2(i)],[y_real1(i) y_real2(i)]);

   plot(x(i),y(i),'-x');
   pause(0);
end