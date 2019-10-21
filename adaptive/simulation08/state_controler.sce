xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t0 = 0;
delta_t = 0.001;
final_time = 10;
t = (0:delta_t:final_time);
t = t'
ifinal=size(t);ifinal=ifinal(1);

F = 10*ones(ifinal, 1);
m = 40;
b = 10;
k = 100;
s = poly(0, 's');
g = syslin('c', 1/m/(s^2 +b*s/m+ k/m))
y = csim(F', t', g)';
filtro = 1/(s+1)^3;
hf = syslin('c', filtro);
u_fil = csim(F', t', hf);
y_fil = csim(y', t', hf);
dy_fil = csim(y', t', s*hf);
d2y_fil = csim(y', t', s*s*hf);

P = 1000000000*eye(3, 3);
theta_plot = zeros(3, ifinal);
p = zeros(ifinal);
m_hat = zeros(ifinal);
k_hat = zeros(ifinal);
b_hat = zeros(ifinal);

p(1) = norm(P, 'fro')
theta = [0;0;0];

entrada = 2*ones(ifinal,1);
x1 = zeros(ifinal);
x2 = zeros(ifinal);
x0=[0;0]; 
A = [-b/m 1; -k/m 0];
A_control = [0 1; -k/m -b/m];
B_control = [0; 1];
C_control = [-1/m 0];
B = [0; 1/m];    
C = [1 0];    
x1(1) = x0(1);
x2(1) = x0(2);
K_control = [10 9.75];
x_alpha = [0.2; 0];
x_control = [0.2; 0];
x1_alpha = zeros(ifinal);
x2_alpha = zeros(ifinal);
x1_alpha(1) = x_alpha(1);
x2_alpha(1) = x_alpha(2);
y_control = zeros(ifinal);
for i=1:ifinal
    fi = [u_fil(i) -y_fil(i) -dy_fil(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(d2y_fil(i) - fi'*theta);
    d2y_fil_hat(i) = fi'*theta;
    theta_plot(:,i) = theta;
    xp = A*x0 + B*(entrada(i) - K_control*x_alpha);
    x0 = x0 + xp*delta_t
    x1(i) = x0(1)
    x2(i) = x0(2)
    A_p = [-theta_plot(3, i) 1; -theta_plot(2, i) 0];
    b_p = [0;theta_plot(1, i)];
    K = [44;10] - [theta_plot(3,i); theta_plot(2,i)];  
    y_alpha = [1 0]*x_alpha; 
    xdot_alpha = A_p*x_alpha + b_p*(entrada(i) - K_control*x_alpha) + K*(x1(i)-y_alpha(1))
    x_alpha_n = x_alpha + xdot_alpha*delta_t;    
    x1_alpha(i) = x_alpha(1);
    x2_alpha(i) = x_alpha(2);
    A_control = [0 1;-theta_plot(2,i)  -theta_plot(3,i)];
    B_control = [0;1];
    xdot_control = A*x_control + B*(entrada(i) - K_control*x_alpha);
    x_alpha = x_alpha_n
    x_control = x_control + xdot_control*delta_t;    
    y_control(i) = C*x_control;   
end


figure(1)
plot(t,theta_plot')


figure(4); //Criando uma janela gráfica 
plot2d(t, x1,5);//vermelho
plot2d(t, x1_alpha)

figure(5); //Criando uma janela gráfica 
plot2d(t, x2',5);//vermelho
plot2d(t, x2_alpha);



figure(6); //Criando uma janela gráfica 
//plot2d(t, entrada, 5);//vermelho
plot2d(t, y_control)


