//Diretório de trabalho
// cd /home/claudiomachado/Documentos/IFSUL/EngEletrica/ControleNaoLinear/Simulacao/Scilab


xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console


//Dinâmica do sistema de segunda ordem
// Equação de estados
//dxdt=[x2;-(K/M)*x1-(B/M)*x2+(1/M)*F];
//x1=x(1,:);//[m]Posição 
//x2=x(2,:);//[m/s]Velocidade

// Parâmetros do modelo:
F=5;//força (N).
M=1;//Massa (kg).
B=0.5;//Coeficiente de atrito viscoso (Ns/m).
K=1;//Rigidês elástica de mola (N/m).

// Matrizes A e b da equação de estado:
////dxdt=A*x+b*F;// Equação de estado
A=[0,1;-K/M,-B/M];
b=[0;1/M];
//
//// Condições iniciais:
t0=0; // Tempo inicial de simulação (s).
x10=[0.0]; // Condição inicial de posição.
x20=[0.0]; // Condição inicial de velocidade.
x0=[x10;x20];

// Parâmetros de simulação:
Ts=1/100;//[s] Período de integração numérica.
tmax=50;//[s] Duração da simulação.
t = [t0:Ts:tmax]; //Vetor de tempo da simulação.

// Definição da função a ser integrada (dx/dt)
function dxdt=f(t,x)
x1=x(1,:)
x2=x(2,:)
dxdt=[x2;-(K/M)*x1-(B/M)*x2+(1/M)*F];// Equação de estados
//dxdt=A*x+b*F;// Equação de estados
endfunction


// Resolve a equação diferencial via Runge-Kutta RKF45:
x = ode("rkf",x0,t0,t,f);
// Desenha o gráfico da solução:
x1=x(1,:);
x2=x(2,:);
plot2d(t,x1);

a=gca();//obtendo as configurações dos eixos
a.parent.background= 8;// background = branco
a.grid = [2,2];//o valor [5,5] dá a cor azul do grid
//a.x_location ="middle"; 
//a.y_location ="left"; 
a.title.text="Simulação no SCILAB"
a.title.font_size = 4;
a.x_label.text="Tempo [s]";
a.x_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
a.x_label.font_size = 4;
a.y_label.text="Posição [m]";
a.y_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
a.y_label.font_size = 4;
a.thickness = 1; //espessura da linha do grid e dos eixos
a.font_size = 3;
//a.data_bounds = [xmin,ymin;xmax,ymax]
