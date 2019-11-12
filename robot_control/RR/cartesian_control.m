clear 
clc
close all
m1 = 1; 
m2 = 1; 
L1 = 2;
L2 = 2;
delta_t = 0.005;
final_time=20;
t = 0:delta_t:final_time;
t = t';
t_plot = 0:delta_t:(final_time-2*delta_t);
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
%theta1_d = 0.1*sin(2*pi*t/2);
%theta2_d = 0.1*cos(2*pi*t/2);

%x = cos(theta1_d) + cos(theta1_d + theta2_d);
%y = sin(theta1_d) + sin(theta1_d + theta2_d);

%y = 1*ones(sizes(1),1);
%x = 1*sin(5*t);

x_d = 0.5*cos(10*t);
y_d = 2*ones(sizes(1),1);
cos_theta2 = (x_d.^2 + y_d.^2 - L1^2*ones(sizes(1),1) - L2^2*ones(sizes(1),1))/(2*L1*L2);
sin_theta2 = (1 - cos_theta2.^2).^0.5;
theta2_d = atan2(sin_theta2, cos_theta2);
theta1_d = atan2(y_d, x_d) - atan2(L2*sin(theta2_d), L1*ones(sizes(1),1) + L2*cos(theta2_d));


x_dotd = (1/delta_t)*diff(x_d);
y_dotd = (1/delta_t)*diff(y_d);
x_dot2d = (1/delta_t)*diff(x_dotd);
y_dot2d = (1/delta_t)*diff(y_dotd);



y_des = [x_d y_d]';
y_dot_des = [x_dotd y_dotd]';
y_dot2_des = [x_dot2d y_dot2d]';

q0 = [0; 0.5];
y0 = [L1*cos(q0(1))+L2*cos(q0(1)+q0(2));L1*sin(q0(1))+L2*sin(q0(1)+q0(2))];
q_dot0 = [0; 0];
gain_p = 1000;
gain_v = 2*gain_p^0.5;
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];

g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);

error_old = (y_des(:,1) - y0);
error_plot = zeros(sizes(1)-2,2);




for i=1:sizes(1)-2
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q0(2)) m2*L2^2+m2*L1*L2*cos(q0(2)); 
         m2*L2^2 + m2*L1*L2*cos(q0(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot0(1)*q_dot0(2)+q_dot0(2)^2)*sin(q0(2));
          m2*L1*L2*q_dot0(1)^2*sin(q0(2))];
    G = [(m1+m2)*g*L1*cos(q0(1))+m2*g*L2*cos(q0(1)+q0(2)); m2*g*L2*cos(q0(1)+q0(2))];   
    J = [-L1*sin(q0(1))-L2*sin(q0(1)+q0(2)) -L2*sin(q0(1)+q0(2));
        L1*cos(q0(1))+L2*cos(q0(1)+q0(2)) L2*cos(q0(1)+q0(2))];
    J_ponto = [-L1*cos(q0(1))-L2*cos(q0(1)+q0(2)) -L2*cos(q0(1)+q0(2));
        -L1*sin(q0(1))-L2*sin(q0(1)+q0(2)) -L2*sin(q0(1)+q0(2))];
    x_real  = L1*cos(q0(1)) + L2*cos(q0(1)+q0(2));
    y_real  = L1*sin(q0(1)) + L2*sin(q0(1)+q0(2));
    error_new = (y_des(:,i) - [x_real;y_real]);
    cart_velocity = J*q_dot0;    
    u = k_v.*(y_dot_des(:,i)-cart_velocity) + k_p.*error_new;
    if rank(J)<2
      tau = [0;0];
    else
      tau = M*inv(J)*(y_dot2_des(:,i) - J_ponto*q_dot0 + u) + V +G;
    endif
    q_dot2 = inv(M)*(tau-V-G);
    q_dot = q_dot0+ delta_t*q_dot2;
    q = q0 + delta_t*q_dot + 0.5*(q_dot2)*(delta_t)^2;
    error_plot(i,:) = error_new;
    q0 = q;
    q_dot0 = q_dot;
    computed_torque1(i) = tau(1);
    computed_torque2(i) = tau(2);
    q_plot(i,:) = q;
end

figure(1)
subplot(2,1,1);
plot(t, x_d);
title('X - Axis')

subplot(2,1,2);
plot(t, y_d);
title('Y - Axis')

figure(2)
subplot(2,1,1);
plot(t, theta1_d);
title('Angle - First Joint')

subplot(2,1,2);
plot(t, theta2_d);
title('Angle -  Second Joint')

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

trajectory_x = [];
trajectory_y = [];

for i=1:sizes(1)-2
  figure(7)

   hold off
   plot([0 x_real1(i)],[0 y_real1(i)]);
   

   axis([-1 3  -1 3]);
   grid on;

   title('RR Robot Arm');
   hold on
   xlabel('x (m)');
   plot(x_d(i),y_d(i),'-x');

   plot([x_real1(i) x_real2(i)],[y_real1(i) y_real2(i)]);
   trajectory_x = [trajectory_x; x_real2(i)];
   trajectory_y = [trajectory_y; y_real2(i)];
   plot(trajectory_x, trajectory_y);
   pause(0);
end