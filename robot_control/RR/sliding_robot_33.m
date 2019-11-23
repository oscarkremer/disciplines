clear 
clc
close all
m1 = 1; 
m2 = 1; 
L1 = 1;
L2 = 1;
delta_t = 0.001;
final_time = 10
t = 0:delta_t:final_time;
t = t';
t_plot = 0:delta_t:final_time-2*delta_t;
t_plot = t_plot';
sizes = size(t);

y = 0.3*sin(5*t) + 0.4*ones(sizes(1),1);
x = 0.3*cos(5*t) + 0.2*ones(sizes(1),1);

cos_theta2 = (x.^2 + y.^2 - L1^2*ones(sizes(1),1) - L2^2*ones(sizes(1),1))/(2*L1*L2);
sin_theta2 = (1-cos_theta2.^2).^0.5;
theta2_d = atan2(sin_theta2, cos_theta2);
theta1_d = atan2(y, x) - atan2(L2*sin(theta2_d), L1*ones(sizes(1),1) + L2*cos(theta2_d));


theta1_dotd = (1/delta_t)*diff(theta1_d);
theta2_dotd = (1/delta_t)*diff(theta2_d);
theta1_dot2d = (1/delta_t)*diff(theta1_dotd);
theta2_dot2d = (1/delta_t)*diff(theta2_dotd);

q_des = [theta1_d theta2_d]';
q_dot_des = [theta1_dotd theta2_dotd]';
q_dot2_des = [theta1_dot2d theta2_dot2d]';

q = [0.1; 0];
q_dot = [0; 0];

g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_torque1 = zeros(sizes(1)-2,2);
computed_torque2 = zeros(sizes(1)-2,2);
lambda = 5*eye(2);
K = 75*eye(2);
sigma = 0.1;
error_plot = zeros(sizes(1)-2,2);
for i=1:sizes(1)-2
    error = (q_des(:,i) - q);
    error_dot = q_dot_des(:,i) - q_dot;
    M = [(m1+m2)*L1^2+m2*L2^2+2*m2*L1*L2*cos(q(2)) m2*L2^2+m2*L1*L2*cos(q(2)); 
         m2*L2^2 + m2*L1*L2*cos(q(2)) m2*L2^2]; 
    V = [-m2*L1*L2*(2*q_dot(1)*q_dot(2) + q_dot(2)^2)*sin(q(2));
          m2*L1*L2*q_dot(1)^2*sin(q(2))];
    Vm = [-2*m2*L1*L2*q_dot(2)*sin(q(2)) -m2*L1*L2*q_dot(2)*sin(q(2));
          m2*L1*L2*q_dot(1)*sin(q(2)) 0];
    G = [(m1+m2)*g*L1*cos(q(1))+m2*g*L2*cos(q(1)+q(2)); m2*g*L2*cos(q(1)+q(2))];   
    r = lambda*error + error_dot;
    if abs(r(1)) > sigma
        r(1) = sign(r(1));       
    else
        r(1) = r(1)/sigma;
    end
    if abs(r(2)) > sigma
        r(2) = sign(r(2));       
    else
        r(2) = r(2)/sigma;
    end
    q_dotr = lambda*error + q_dot;
    q_2dotr = lambda*error_dot + q_dot2_des(:,i);
    tau = M*q_2dotr + Vm*q_dotr + G + K*r;
    q_dot2 = inv(M)*(tau-V-G);
    q_dot = q_dot+ delta_t*q_dot2;
    q = q + delta_t*q_dot + 0.5*(q_dot2)*(delta_t)^2;
    error_plot(i,:) = error;
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


for i=1:sizes(1)-2
   figure(7)
   plot(x(i),y(i),'-x');

   hold off
   plot([0 x_real1(i)],[0 y_real1(i)]);

   axis([-2 2  -0.5 2.0]);
   grid on;

   title('RP Robot Arm');
   hold on
   xlabel('x (m)');
   plot([x_real1(i) x_real2(i)],[y_real1(i) y_real2(i)]);

   pause(0);
end