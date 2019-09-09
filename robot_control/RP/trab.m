clear 
clc
close all
m = 1;
t = 0:0.0005:6;
t = t';
t_plot = 0:0.0005:5.9990;
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
%y = sin(5*t) + 1.1*ones(sizes(1),1);
%x = cos(5*t);

y = ones(sizes(1),1);
x = sin(5*t);
r_d = (x.^2 + y.^2).^0.5;
theta_d = atan2(y, x);
r_dotd = diff(r_d);
theta_dotd = diff(theta_d);
r_dot2d = diff(r_dotd);
theta_dot2d = diff(theta_dotd);

q_des = [theta_d r_d]';
q_dot_des = [theta_dotd r_dotd]';
q_dot2_des = [theta_dot2d r_dot2d]';

q0 = [3.141592/2; 1.5];
q_dot0 = [0; 0];
gain_p = 5000;
gain_v = 2*(gain_p^0.5);
k_p = [gain_p; gain_p];
k_v = [gain_v; gain_v];
g = 9.8;
q_plot = zeros(sizes(1)-2,2);
computed_force = zeros(sizes(1)-2,2);
computed_torque = zeros(sizes(1)-2,2);

error_plot = zeros(sizes(1)-2,2);

tau = [0;0];
for i=1:sizes(1)-2
    M = [m*q0(2)*q0(2) 0; 
         0 m]; 
    V = [2*m*q0(2)*q_dot0(2)*q_dot0(1);
          -m*q0(2)*(q_dot0(1)^2)];
    G = [m*g*q0(2)*cos(q0(1)); m*g*sin(q0(1))];   
    error_new = (q_des(:,i) - q0);
    q_dot2 = q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0);
    q_dot = q_dot0+ 0.0005*q_dot2;
    q = q0 + 0.0005*q_dot + 0.5*(q_dot2)*(0.0005)^2;
    q0 = q;
    q_dot0 = q_dot;
    tau = M*(q_dot2_des(:,i) + k_v.*(q_dot_des(:,i)-q_dot0) + k_p.*(q_des(:,i) - q0)) + V +G;
    computed_torque(i) = tau(1);
    computed_force(i) = tau(2);
    error_plot(i,:) = error_new;
    q_plot(i,:) = q;
end

x_real  = q_plot(:,2).*cos(q_plot(:,1));
y_real  = q_plot(:,2).*sin(q_plot(:,1));


for i=1:sizes(1)-2
   figure(1)
   hold off
   plot([0 x_real(i)],[0 y_real(i)])
   axis([-2.2 2.2  0 1.5]);
   grid on;
   title('RP Robot Arm');
   hold on
   xlabel('x (m)');
   plot(x(i),y(i),'-x');
   pause(0);
end