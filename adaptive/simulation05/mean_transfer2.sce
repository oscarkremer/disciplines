clear
close

t = (0:0.05:10);
t = t'
ifinal=size(t);ifinal=ifinal(1);
F = sin(4*t);
m = 40;
b = 10;
k = 100;
s = poly(0, 's');
g = syslin('c', 1/(m*s^2 +b*s+ k))
y = csim(F', t', g)';
filtro = 1/(s+1)^3;
hf = syslin('c', filtro);
u_fil = csim(F', t', hf);
y_fil = csim(y', t', hf);
dy_fil = csim(y', t', s*hf);
d2y_fil = csim(y', t', s*s*hf);

P = 10000000*eye(3, 3);
theta_plot = zeros(3, ifinal);
p = zeros(ifinal);
p(1) = norm(P, 'fro')
theta = [0;0;0];

for i=2:ifinal
    fi = [u_fil(i) -y_fil(i) -dy_fil(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(d2y_fil(i) - fi'*theta);
    theta_plot(:,i) = theta
end


figure(1)
plot(t,theta_plot')


figure(2)
plot(t,p)
