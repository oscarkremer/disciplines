clear
close
delta_t=0.01
t = (0:delta_t:100);
t = t'
ifinal=size(t);ifinal=ifinal(1);
uc = 0.05*sin(16.1*t);
s = poly(0, 's');
g = syslin('c', 2/(s + 2))
ym = csim(uc', t', g)';
y0 = [0];
y_plot = zeros(ifinal,1)
y_plot(1) = y0;
a = 1;
b = 0.5;
theta1 = 0;
theta2 = 0;
theta_plot1 = zeros(ifinal,1)
theta_plot2 = zeros(ifinal,1)
theta_plot1(1) = theta1;
theta_plot2(1) = theta2;
gamma1 = 200.539
gamma2 = 200.539
error_plot = zeros(ifinal, 1)
ITAE = 0;
model_error = y0-ym(1)
for i=2:ifinal
    theta1 = theta1 + delta_t*(-gamma1*model_error*uc(i-1));
    theta2 = theta2 + delta_t*(gamma2*model_error*y0);
    u = theta1*uc(i-1)-theta2*y0
    y0 = (1-a*delta_t)*y0 + b*delta_t*u
    model_error = (y0 - ym(i))
    ITAE = ITAE + abs(model_error)*t(i)*delta_t
    theta_plot1(i) = theta1;
    theta_plot2(i) = theta2;
    y_plot(i) = y0    
end

figure(1)
plot(t,y_plot,'-k')
plot(t,ym,'-b')


figure(2)
plot(t,theta_plot1)


figure(3)
plot(t,theta_plot2)


figure(4)
plot(t, y_plot - ym)
h = gca(); // get current axes

h.data_bounds = [0, 100 ; 0, 0.6];
