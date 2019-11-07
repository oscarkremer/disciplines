xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t0 = 0;
delta_t = 0.001;
final_time = 10;
t = (0:delta_t:final_time);
t = t'
ifinal=size(t);ifinal=ifinal(1);


m = 40;
b = 10;
k = 100;
P = 1000000*eye(3, 3);
theta_plot = zeros(3, ifinal);
m_hat = zeros(ifinal);
k_hat = zeros(ifinal);
b_hat = zeros(ifinal);
theta = [0;0;0];
xd = sin(2*t);
xpd = 2*cos(2*t);
xp2d = -4*sin(2*t);
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
x_alpha = [0.2; 0];
xdot_alpha = [0;0];
x1_alpha = zeros(ifinal);
x2_alpha = zeros(ifinal);
x1_alpha(1) = x_alpha(1);
x2_alpha(1) = x_alpha(2);
xp = [x0(1);0];
kp = 100;
kd = 20;

for i=1:ifinal
    control =  m*(xp2d(i)+kp*(xd(i)-x_alpha(1)) + kd*(xpd(i)-xdot_alpha(1))) + k*x_alpha(1)+b*xdot_alpha(1)
    fi = [control -x0(1) -xp(1)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    theta = theta + K*((1/m)*(control -b*xp(1)-k*x0(1)) - fi'*theta);
    theta_plot(:,i) = theta;
    A_p = [-theta_plot(3, i) 1; -theta_plot(2, i) 0];
    b_p = [0;theta_plot(1, i)];
    xp = A*x0 + B*control;
    x0 = x0 + xp*delta_t
    x1(i) = x0(1)
    x2(i) = x0(2)
    K_obs = [44;10] - [theta_plot(3,i); theta_plot(2,i)];  
    y_alpha = [1 0]*x_alpha; 
    xdot_alpha = A_p*x_alpha + b_p*(control) + K_obs*(x1(i)-y_alpha(1))
    x_alpha = x_alpha + xdot_alpha*delta_t;    
    x1_alpha(i) = x_alpha(1);
    x2_alpha(i) = x_alpha(2);
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
plot2d(t, x1)


