clear
close
delta_t=0.001
t = (0:delta_t:10);
t = t'
ifinal=size(t);ifinal=ifinal(1);
F = t;
a = 1;
b = 2;
s = poly(0, 's');
g = syslin('c', a/(s + b))
y = csim(F', t', g)';
filtro = 1/(s+1)^2;
hf = syslin('c', filtro);
u_fil = csim(F', t', hf);
y_fil = csim(y', t', hf);
dy_fil = csim(y', t', s*hf);

P = 10000*eye(2, 2);
theta_plot = zeros(2, ifinal);
p = zeros(ifinal);
p(1) = norm(P, 'fro')
theta = [10;10];
y0 = [0];
right_y0 = y0;
omega_n = 5;
zeta = 0.9;
signal = [ones(2*ifinal/10,1); zeros(2*ifinal/10,1); ones(2*ifinal/10,1); zeros(2*ifinal/10,1); ones(2*ifinal/10 + 1,1)];
y_plot = zeros(ifinal,1)
y_plot(1) = y0;
right_ki = omega_n^2/a
right_kp = (2*zeta*omega_n - b)/a; 
integral_error = 0;
right_integral_error = 0;
for i=2:ifinal
    fi = [u_fil(i) -y_fil(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(2,2) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(dy_fil(i) - fi'*theta);
    theta_plot(:,i) = theta
    ki =  omega_n^2/theta(1)
    kp = (2*zeta*omega_n - theta(2))/theta(1);
    signal_error = signal(i)-y0;
    integral_error = integral_error + delta_t*signal_error;
    uc = ki*integral_error + kp*signal_error
    y0 = (delta_t*a*uc+y0)/(1+b*delta_t);
    right_signal_error = signal(i)-right_y0;
    right_integral_error = right_integral_error + delta_t*right_signal_error;
    uc = right_ki*right_integral_error + right_kp*right_signal_error
    right_y0 = (delta_t*a*uc+right_y0)/(1+b*delta_t);
    y_plot(i) = y0    
    y_right_plot(i) = right_y0
end


figure(1)
plot(t,theta_plot')


figure(2)
plot(t,p)


figure(3)
plot(t,signal)
plot(t,y_right_plot, '-r')
plot(t,y_plot, '-k')
