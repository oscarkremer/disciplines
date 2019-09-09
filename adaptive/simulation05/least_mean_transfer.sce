t = (0:0.01:0.5)';
noisegen(0.01, 5, 0.1);
F = feval(t, Noise);


u1 = sin(t);
u2 = cos(t);
u3 = 3*t.^2;

y = 3*u1+ 2*u2+ 6*u3;

P = 100000*eye(3, 3);
PHI = [sin(t);cos(t);3*t.^2] ;
ifinal=size(t);ifinal=ifinal(1);
theta_plot = zeros(3, ifinal);
p = zeros(ifinal);
p(1) = norm(P, 'fro')
theta = [0;0;0];
for i=2:ifinal
    fi = [u1(i) u2(i) u3(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(y(i) - fi'*theta);
    theta_plot(:,i) = theta
end

figure(0)
plot(t,theta_plot)


figure(1)
plot(t,p)
