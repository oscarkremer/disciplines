clear
close

t = (0:0.001:5)';
ifinal=size(t);ifinal=ifinal(1);
noisegen(0.001, 5, 20);
F = feval(t, Noise);
x0 = 0;
x0_p = 0;
x_f = zeros(ifinal,3)
x_2p = zeros(ifinal);
x_p = zeros(ifinal);
x = zeros(ifinal);
m = 40;
b = 10;
k = 100;
x(1) = x0;
x_p(1) = x0_p;

for i=2:ifinal
    x2p = (1/m)*(F(i)-b*x0_p-k*x0);
    x1p = x0_p + x2p*0.001;
    x(i) = x0 + x0_p*0.001 + (0.5)*x2p*0.001*0.001;    
    x_p(i) = x1p;
    x_2p(i) = x2p;
    x0_p = x1p;
    x0 = x(i);
end

P = 100000*eye(3, 3);
theta_plot = zeros(3, ifinal);
p = zeros(ifinal);
p(1) = norm(P, 'fro')
theta = [0;0;0];
for i=2:ifinal
    fi = [F(i) -x(i) -x_p(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(x_2p(i) - fi'*theta);
    theta_plot(:,i) = theta
end


figure(0)
plot(t,F)


figure(1)
plot(t,theta_plot)


figure(2)
plot(t,p)

