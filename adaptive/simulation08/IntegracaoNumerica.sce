xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console

t0=0; // Tempo inicial de simulação (s).
Ts=0.01;//[s] Período de integração numérica.
tmax=50;//[s] Duração da simulação.
t = [t0:Ts:tmax]; //Vetor de tempo da simulação.
ifinal = size(t); ifinal = ifinal(2);

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

//
////Integração numérica de sistema de segunda ordem
////Dinâmica do sistema de segunda ordem
//// Equação de estados
////dx1dt=;
////dx2dt=;
////dxdt=[dx1dt;dx2dt];
////
////////////////////////////////////////////
////////Parâmetros
amp=20;//amplitude do seno
T=5;//períodp do seno [s]
fr=1/T;//frequência [5Hz]
w=2*%pi*fr;//frequência [rad/s]
kp = 100;
kd = 25;
ki = 40;
// Condições iniciais:
x10=[0.0]; // Condição inicial
x20=[0.0]; // Condição inicial
x30=[1.0]; // Condição inicial
x40=[0.0]; // Condição inicial
x50 = [1-x10];
x0=[x10;x20;x30;x40;x50]; //vetor de condições iniciais

k = 100;
b = 10;
m = 40;
// Definição da função a ser integrada (dx/dt)
theta1 = 1/m;
theta2 = k/m;
theta3 = b/m;
theta_p1 = [0;0];
theta_p2 = [0;0;0;0];
function dxdt=f(t,x)
x1=x(1,:); //renomeando o estado
x2=x(2,:);
x3=x(3,:);
x4=x(4,:);
x5=x(5,:);
d=1;
//theta2 = theta_p2(1);
//theta1 = theta_p2(4);
//print(%io(2), theta1)
//renomeando o estado
fi1 = [-x1(1) x2(1)]';
K1 =P1*fi1/(1+fi1'*P1*fi1);
P1 = (eye(2,2) - K1*fi1')*P1;
theta_p1 = theta_p1 + K1*((-b/m)*x1(1)+x2(1) - fi1'*theta_p1);
theta3 = theta_p1(1);
print(%io(2), theta_p1);
//d=0
//B=[0;d/a]
e = d - x3;
uc = kp*e + kd*(-44*x1 - (-theta3-44)*x3-x4) + ki*x5;
dx1dt=-b/m*x1 + x2;
dx2dt=-k/m*x1 + (1/m)*uc;
dx3dt= 44*x1 + (-theta3-44)*x3 + x4;
dx4dt= 10*x1 + (-theta2-10)*x3 + theta1*uc;
dxdt=[dx1dt;dx2dt;dx3dt;dx4dt;e];
endfunction

//// Resolve a equação diferencial via Runge-Kutta RKF45:
x = ode("rkf",x0,t0,t,f);

//// Desenha o gráfico da solução:
x1=x(1,:); //renomeando o estado
x2=x(2,:); //renomeando o estado
x3=x(3,:); //renomeando o estado
x4=x(4,:); //renomeando o estado

figure(0); //Criando uma janela gráfica 
//plot2d(t,x1,5);//vermelho
plot2d(t, x1,5);//vermelho
plot2d(t, x3);//vermelho


figure(1); //Criando uma janela gráfica 
//plot2d(t,x1,5);//vermelho
plot2d(t, x2,5);//vermelho
plot2d(t, x4);//vermelho

//a=gca();//obtendo as configurações dos eixos
//a.parent.background=8; // background = branco
//a.grid = [2,2];//o valor [5,5] dá a cor azul do grid
////a.isoview='on'; //escalas iguais nos dois eixos
////a.x_location ="middle"; 
////a.y_location ="left"; 
//a.title.text="Integração Numérica"
//a.title.font_size = 4;
//a.x_label.text="Tempo [s]";
//a.x_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//a.x_label.font_size = 4;
//a.y_label.text="x1";
//a.y_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//a.y_label.font_size = 4;
//a.thickness = 1; //espessura da linha do grid e dos eixos
//a.font_size = 3;
//ymin=2*min(x2);ymax=2*max(x2);
//a.data_bounds = [0,ymin;tmax,ymax];//a.data_bounds = [xmin,ymin;xmax,ymax]
//
//figure(1); //Criando uma janela gráfica 
//plot2d(t,x2,2);//azul
//
//a=gca();//obtendo as configurações dos eixos
//a.parent.background=8; // background = branco
//a.grid = [2,2];//o valor [5,5] dá a cor azul do grid
////a.isoview='on'; //escalas iguais nos dois eixos
////a.x_location ="middle"; 
////a.y_location ="left"; 
//a.title.text="Integração Numérica"
//a.title.font_size = 4;
//a.x_label.text="Tempo [s]";
//a.x_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//a.x_label.font_size = 4;
//a.y_label.text="x2";
//a.y_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//a.y_label.font_size = 4;
//a.thickness = 1; //espessura da linha do grid e dos eixos
//a.font_size = 3;
//ymin=2*min(x2);ymax=2*max(x2);
//a.data_bounds = [0,ymin;tmax,ymax];//a.data_bounds = [xmin,ymin;xmax,ymax]
//
