//Diretório de trabalho
// cd /home/claudiomachado/Documentos/IFSUL/EngEletrica/ControleNaoLinear/Simulacao/Scilab

clc();
xdel(winsid());
clear();
//stacksize('max');

// Parâmetros do modelo:
J=1.6/1000/100/100; //Momento de inércia do motor [gcm^2]
R=5.8;//Resistencia elétrica do motor [ohm].
L=135D-6;//Indutância do motor [uH].
Ki=1.496D-3;//constante de força contra-eletromotriz do motor, dada em [mV/min^-1] = [mV/RPM]
Kt=14.29D-3;//constante de torque, dada em [mNm/A]
Tf=0.4D-3;// Torque de atrito do motor [mNm]

// Condições iniciais:
x1d = 700;
x01=[0.0];//
x02=[0.0];//
x03=[x1d-x01];//

// Parâmetros de simulação:
t0=0; // Tempo inicial (s).
Ts=1/10000;// Período (s).
tmax=1;//5*Ts;// Duração da simulação (s).
t = [t0:Ts:tmax]; // Gera o tempo.
// Definição da função a ser integrada (dx/dt)
kp = 0.001;
kd =0.01;
ki = 0.010;
function dxdt=f(t,x)
    x1=x(1,:);
    x2=x(2,:);
    x3=x(3,:);
    if abs(x2)<0.01 then
        if abs(Kt*x2)<Tf then
        Bf=Kt*x2;
        else
        Bf=Tf;
        end      
    else        
        Bf=Tf;
    end
    e = x1d - x1;
    u = kp*e + ki*x3;
    dxdt=[-(Bf/J)+(Kt*x2/J);-(R*x2/L)-(Ki*x1/L) + u/L; e];// Equação de estado

endfunction

// Resolve a equação diferencial via Runge-Kutta RKF45:
x = ode("rkf",[x01;x02;x03],t0,t,f);
    x1=x(1,:);
    x2=x(2,:);
    x3=x(3,:);
//// Desenha o gráfico da solução:
figure(0)
plot2d(t,x1);

figure(1)
plot2d(t,x2);

figure(2)
plot2d(t,x3);


