clear
close
delta_t=0.001
t = (0:delta_t:20);
t = t'
ifinal=size(t);ifinal=ifinal(1);
uc = 1.5*cos(t);
s = poly(0, 's');
g = syslin('c', 2/(s + 1))
ym = csim(uc', t', g)';
y0 = [0];
y_plot = zeros(ifinal,1)
y_plot(1) = y0;
theta = 0;
theta_plot = zeros(ifinal,1)
theta_plot(1) = theta;
gamma_m = 1;
for i=2:ifinal
    model_error = (y0 - ym(i))
    theta = theta + delta_t*(-gamma_m*model_error*ym(i));
    theta_plot(i) = theta;
    u = theta*uc(i)
    y0 = (1-delta_t)*y0 + delta_t*u
    y_plot(i) = y0    
end

figure(1)
plot(t,y_plot,'-k')
plot(t,ym,'-b')


figure(2)
plot(t,theta_plot)


