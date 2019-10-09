//Diretório de trabalho
// cd /home/claudiomachado/Documentos/IFSUL/EngEletrica/ControleNaoLinear/Simulacao/Scilab

clc();
xdel(winsid());
clear();
//stacksize('max');

// Parâmetros do modelo:
J=1.8/1000/100/100; //Momento de inércia do motor [gcm^2]
R=5.8;//Resistencia elétrica do motor [ohm].
L=135D-6;//Indutância do motor [uH].
Ki=1.496D-3;//constante de força contra-eletromotriz do motor, dada em [mV/min^-1] = [mV/RPM]
Kt=14.29D-3;//constante de torque, dada em [mNm/A]
Tf=0.4D-3;// Torque de atrito do motor [mNm]
B=Kt/J;
C=Ki/L;
D=R/L;
E=1/L;

//controlador
w=500;//velocidade desejada [rpm]
ow=600;//offset de velocidade [rpm]
Tw=2;//periodo de oscilação da velocidade desejada [s]
fw=1/Tw;//frequência de oscilação da velocidade desejada [Hz]
wf=2*%pi*fw;//frequência de oscilação da velocidade desejada [rad/s]
kp=1;//ganho do controlador proporcional


// Condições iniciais:
x01=[0.0];//
x02=[0.0];//

// Parâmetros de simulação:
t0=0; // Tempo inicial (s).
Ts=1/10000;// Período (s).
tmax=5;//5*Ts;// Duração da simulação (s).
t = [t0:Ts:tmax]; // Gera o tempo.
// Definição da função a ser integrada (dx/dt)

function dxdt=f(t,x)
    x1=x(1,:);
    x2=x(2,:);
    if abs(x2)<0.01 then
        if abs(Kt*x2)<Tf then
        Bf=Kt*x2;
        else
        Bf=Tf;
        end      
    else        
        Bf=Tf;
    end
    A=Bf/J;
//    u=5;
    x1d=6;
    e=x1-x1d;
    
    u=-kp*e-1*(-A+B*x2);
    dxdt=[-(A)+(B*x2);-(C*x1)-(D*x2)+(E*u)];// Equação de estado
endfunction

// Resolve a equação diferencial via Runge-Kutta RKF45:
x = ode("rkf",[x01;x02],t0,t,f);
    x1=x(1,:);
    x2=x(2,:);
//// Desenha o gráfico da solução:
figure(0)
plot2d(t,x1);

figure(1)
plot2d(t,x2);


