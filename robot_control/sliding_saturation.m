clear 
clc
close all


delta_t = 0.0002;
sampling_rate = 0.001;
sampling_int = sampling_rate/delta_t;
final_t = 4;
t = 0:delta_t:final_t;
t = t';
t_plot = 0:delta_t:final_t;
t_plot = t_plot';
sizes = size(t);

xd = sin(pi*t/2);
xd_dot = (pi/2)*cos(pi*t/2);
xd_2dot = -(pi/2)^2*sin(pi*t/2);

x = 0;
x_dot = pi/2;

lambda = 20;
eta = 0.1;
phi = 0.1;
x_plot = zeros(sizes(1),1);
u_plot = zeros(sizes(1),1);
error_plot = zeros(sizes(1),1);

for i=1:sizes(1)
    x_plot(i) = x;
    a = abs(sin(t(i))) + 1;
    error = x - xd(i);
    error_dot = x_dot - xd_dot(i);    
    k = 0.5*x_dot^2*abs(cos(3*x))+eta;
    s = error_dot + lambda*error;
    error_plot(i) = error;
    if rem(i-1, sampling_int) == 0
      f_hat = -1.5*x_dot^2*cos(3*x);
      u_hat = -f_hat + xd_2dot(i) - lambda*error_dot;
      if abs(s)>phi
        u = u_hat - k*sign(s/phoi);
      else
        u = u_hat - k*s/phi;        
      end  
    else
      u  = u;
    end
    u_plot(i) = u;
    x_2dot = -a*x_dot^2*cos(3*x) + u;
    x_dot = x_dot + delta_t*x_2dot;
    x = x + x_dot*delta_t + 0.5*x_2dot*(delta_t)^2;
end



figure(1)
plot(t, x_plot, t, xd);
title('Trajectory Following')



figure(4)
plot(t, u_plot);
title('Control Signal')


figure(6)
plot(t, error_plot);
title('Tracking Error')
