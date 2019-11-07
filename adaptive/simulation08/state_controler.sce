
xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t0=0; // Tempo inicial de simulação (s).
Ts=1/250;//[s] Período de integração numérica.
tmax=50;
t = [t0:Ts:tmax]; //Vetor de tempo da simulação.

ifinal=size(t');ifinal=ifinal(1);

F = ones(ifinal, 1);

m = 40;
b = 10;
k = 100;
s = poly(0, 's');
g = syslin('c', 1/(m*s^2 +b*s+ k))
y = csim(F', t, g)';
filtro = 1/(s+2)^3;
hf = syslin('c', filtro);
u_fil = csim(F', t', hf);
y_fil = csim(y', t', hf);
dy_fil = csim(y', t', s*hf);
d2y_fil = csim(y', t', s*s*hf);

P = 100000*eye(3, 3);
p = zeros(ifinal);
p(1) = norm(P, 'fro')
theta = [0;0;0];

for i=2:ifinal
    fi = [u_fil(i) -y_fil(i) -dy_fil(i)]';
    K =P*fi/(1+fi'*P*fi);
    P = (eye(3,3) - K*fi')*P;
    p(i) = norm(P, 'fro')
    theta = theta + K*(d2y_fil(i) - fi'*theta);
end

theta_p1 = [0;0;0;0];
theta_p2 = [0;0];
theta1 = 1/m;
theta2 = k/m;
theta3 = b/m;
k1 =44
function dxdt=f(t,x)
    x1=x(1,:); //renomeando o estado
    x2=x(2,:); //renomeando o estado
    x3=x(3,:);
    x4=x(4,:);
    d=5;
    dx1dt=-b/m*x1 + x2;
    dx2dt=-k/m*x1 + (50/m)*x3 + (10/m)*x4 + (1/m)*d;
    dx3dt= k1*x1 + (-theta3-k1)*x3 + x4;
    dx4dt= 10*x1 + (-theta2-10+50*theta1)*x3 + 10*theta1*x4 + theta1*d;
    dxdt=[dx1dt; dx2dt; dx3dt; dx4dt];
endfunction

x10=[1.0]; // Condição inicial
x20=[0.0]; // Condição inicial
x30=[1.0]; // Condição inicial
x40=[0.0]; // Condição inicial
x0=[x10;x20;x30;x40]

x = ode("rkf",x0,t0,t,f);

//// Desenha o gráfico da solução:
x1=x(1,:); //renomeando o estado
x2=x(2,:); //renomeando o estado
figure(0); //Criando uma janela gráfica 
//plot2d(t,x1,5);//vermelho
plot2d(t, x1,5);//vermelho
