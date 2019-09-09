//Diretório de trabalho
// cd /home/claudio/Dropbox/IFSUL/EngenhariaEletrica/ControleNaoLinear/Simulacao


//Objetivo
//1. Uso do Scilab
//2. Plotagens
//3. Funções
//4. Vetores
//5. Integração numérica

//////////////////////////////////////
xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
//////////////////////////////////////

// Vetor de tempo
t0=0; // Tempo inicial de simulação (s).
Ts=1/250;//[s] Período de integração numérica.
tmax=50;//[s] Duração da simulação.
t = [t0:Ts:tmax]; //Vetor de tempo da simulação.

////////////////////////////////////////////
////////Plotagem
//amp=20;//amplitude do seno
//T=5;//períodp do seno [s]
//fr=1/T;//frequência [Hz]
//w=2*%pi*fr;//frequência [rad/s]
//ys=amp*sin(w*t)
//yc=amp*cos(w*t)
//plot2d(t,ys,5) // Função plotada em vermelho
//plot2d(t,yc,13) // Função plotada em verde
//
////getcolor(); // Padrão de cores de plotagem
////
//
///////////////////////////////////////////////
////// Formatando o gráfico
//a=gca();//obtendo as configurações dos eixos
//a.parent.background=8; // background = branco
//a.grid = [2,2];//o valor [5,5] dá a cor vermelha no grid
////a.isoview='on'; //escalas iguais nos dois eixos
////a.x_location ="middle"; 
////a.y_location ="left"; 
//a.title.text="Seno e cosseno"
//a.title.font_size = 4;
//a.x_label.text="Tempo [s]";
//a.x_label.font_style = 1; //fonte // O TIPO DE FONTE DEPENDE DA CODIFICAÇÃO DE TEXTO UTILIZADA
//a.x_label.font_size = 4;
//a.y_label.text="ys,yc";
//a.y_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//a.y_label.font_size = 4;
//a.thickness = 1; //espessura da linha do grid e dos eixos
//a.font_size = 3;
//ymin=2*min(x2);ymax=2*max(x2);
//a.data_bounds = [0,ymin;tmax,ymax];//a.data_bounds = [xmin,ymin;xmax,ymax]
//

////////////////////////////////////////
////Avaliaçao da derivada de uma função no ponto considerado (ti,yi)
//amp=20;//amplitude do seno
//T=5;//períodp do seno [s]
//fr=1/T;//frequência [Hz]
//w=2*%pi*fr;//frequência [rad/s]
//ys=amp*sin(w*t)
//plot2d(t,ys,5) // Função plotada em vermelho
//ti=T/3 //ponto considerado na derivada
//yi=amp*sin(w*ti)
//plot2d(ti,yi,-9)
////derivada da função seno
//dx2dt=amp*w*cos(w*ti)
////reta tangente ao seno no ponto considerado (ti,yi)
//y=dx2dt*(t-ti)+yi
//plot2d(t,y,1)
//
////Formatação do gráfico
//a=gca();//obtendo as configurações dos eixos
//a.parent.background=8; // background = branco
//a.grid = [2,2];//o valor [5,5] dá a cor azul do grid
//a.x_location ="middle"; //posição do eixo x
//a.y_location ="left";  //posição do eixo x
//ymin=1.5*min(ys);ymax=1.5*max(ys);
//a.data_bounds = [0,ymin;tmax,ymax];//a.data_bounds = [xmin,ymin;xmax,ymax]
//
////Integração numérica
//// Equação de estados
////dx1dt=;
////dx2dt=;
////dx3dt=;
////dxdt=[dx1dt;dx2dt;dx3dt];
//
////////////////////////////////////////////
////////Parâmetros
amp=20;//amplitude do seno
T=5;//períodp do seno [s]
fr=1/T;//frequência [Hz]
w=2*%pi*fr;//frequência [rad/s]


// Condições iniciais:
x10=[-5.0]; // Condição inicial
x20=[10.0]; // Condição inicial
x30=[30.0]; // Condição inicial
x40=[10.0]; // Condição inicial
x0=[x10;x20;x30;x40]; //vetor de condições iniciais

////
//////// Definição da função a ser integrada (dx/dt)
//function dxdt=f(t,x)
//x1=x(1,:); //renomeando o estado
//x2=x(2,:); //renomeando o estado
//x3=x(3,:); //renomeando o estado
//x4=x(4,:); //renomeando o estado
//dx1dt=-0.2*x1+0.3;
//dx2dt=amp*w*cos(w*t);
//dx3dt=-0.5;
//dx4dt=0;
//dxdt=[dx1dt;dx2dt;dx3dt;dx4dt];// Equação de estados
//endfunction
//
//// Resolve a equação diferencial via Runge-Kutta RKF45:
//x = ode("rkf",x0,t0,t,f);
//
////// Desenha o gráfico da solução:
//x1=x(1,:); //renomeando o estado
//x2=x(2,:); //renomeando o estado
//x3=x(3,:); //renomeando o estado
//x4=x(4,:); //renomeando o estado
//figure(0); //Criando uma janela gráfica 
//plot2d(t,x1,5);//vermelho
//plot2d(t,x2,2);//azul
//plot2d(t,x3,25);//vermelho
//plot2d(t,x4,13);//verde
//
//////////Analisando o resultado 
////// seno
//plot2d(t,x20+a*sin((2*%pi/T)*t),1)
//////ponto considerado
//ti=5.5
//yi=x20+a*sin(w*ti)
//plot2d(ti,yi,-9)
////valor da derivada da função seno no ponto
//dx2dt=20*w*cos(w*ti)
////reta tangente ao seno
//y=dx2dt*(t-ti)+yi
//plot2d(t,y,1)
//////
//////a=gca();//obtendo as configurações dos eixos
//////a.parent.background=8; // background = branco
//////a.grid = [2,2];//o valor [5,5] dá a cor azul do grid
////////a.isoview='on'; //escalas iguais nos dois eixos
////////a.x_location ="middle"; 
////////a.y_location ="left"; 
//////a.title.text="Integração Numérica"
//////a.title.font_size = 4;
//////a.x_label.text="Tempo [s]";
//////a.x_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//////a.x_label.font_size = 4;
//////a.y_label.text="x1,x2,x3,x4";
//////a.y_label.font_style = 1; //tipo de fonte - 1=normal; 3=itálico
//////a.y_label.font_size = 4;
//////a.thickness = 1; //espessura da linha do grid e dos eixos
//////a.font_size = 3;
//////ymin=2*min(x2);ymax=2*max(x2);
//////a.data_bounds = [0,ymin;tmax,ymax];//a.data_bounds = [xmin,ymin;xmax,ymax]
//
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
a=1;
b=1;
c=1;
A=[0,1;-c/a,-b/a]


// Condições iniciais:
x10=[1.0]; // Condição inicial
x20=[0.0]; // Condição inicial
x0=[x10;x20]; //vetor de condições iniciais


// Definição da função a ser integrada (dx/dt)
function dxdt=f(t,x)
x1=x(1,:); //renomeando o estado
x2=x(2,:); //renomeando o estado
//d=amp*w*cos(w*t);
d=0
//B=[0;d/a]
dx1dt=x2;
dx2dt=(-c*x1-b*x2+d)/a;
dxdt=[dx1dt;dx2dt];
//dxdt=A*x+B
endfunction

//// Resolve a equação diferencial via Runge-Kutta RKF45:
x = ode("rkf",x0,t0,t,f);

//// Desenha o gráfico da solução:
x1=x(1,:); //renomeando o estado
x2=x(2,:); //renomeando o estado
figure(0); //Criando uma janela gráfica 
//plot2d(t,x1,5);//vermelho
plot2d(x1,x2,5);//vermelho

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
