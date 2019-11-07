%% Rotina para implementacao do algoritmo 3.2 NO EXEMPLO 3.5
clear all, clc, close all

%% TEMPO
dt=1e-3; TMAX = 100;
t=dt:dt:TMAX;
N = length(t);

%% Dados originais
%A = [1 -1.6065 0.6065];
%B = [0.1065 0.0902];
% Estes dados ser�o agora estimados

%% Especifica��es
Am = [1 -1.3205 0.4966];
Bm = [0.1761 0];
A0 = [1 0];

% Condi��es de compatibilidade
%length(A) = length(Am)
%length(B) = length(Bm)
%length(A0) = length(A)-length(Bplus)-1
%length(A0) = d0-1
%Bm = conv(Bminus,Bml)

%[R1,S1,T1] = mdpp(A,B,Am,Bm,A0,'c')

%[R2,S2,T2] = mdpp(A,B,Am,Bm,A0,'s')

%% C�LCULO RECURSIVO
% Vari�veis
%P = 1e7*eye(4);
P = [100 0 0 0; 0 100 0 0; 0 0 1 0; 0 0 0 1];
THETA=[0;0;0.01;0.2];
a1(1)=THETA(1,1);  
a1(2)=THETA(1,1); 
a2(1)=THETA(2,1);  
a2(2)=THETA(2,1); 
b0(1)=THETA(3,1);  
b0(2)=THETA(3,1); 
b1(1)=THETA(4,1); 
b1(2)=THETA(4,1); 
y(1) = 0;
y(2) = 0;
uc(1) = 0;
uc(2) = 0;
u(1) = 0;
u(2) = 0;

for i=3:(TMAX)
    % Entrada 'uc'
    if i<=25
        uc(i) = 1;
    end
    if i>25 && i<=50
        uc(i) = -1;
    end
    if i>50 && i<=75
        uc(i) = 1;
    end 
    if i>75
        uc(i) = -1;
    end
       
    % Sa�da
    y(i) = -(-1.6065)*y(i-1)-(0.6065)*y(i-2)+(0.1065)*u(i-1)+(0.0902)*u(i-2);
    
    % RLS
    phi = [-y(i-1);-y(i-2);u(i-1);u(i-2)];
    phiT=phi';
    K=P*phi*inv(1+phiT*P*phi);
    P=((eye(4)-K*phiT)*P);
    THETA=THETA+K*(y(i)-phiT*THETA);
    a1(i)=THETA(1,1);
    a2(i)=THETA(2,1);
    b0(i)=THETA(3,1); 
    b1(i)=THETA(4,1); 
    
    A = [1 a1(i) a2(i)];
    B = [b0(i) b1(i)];
    [R,S,T] = mdpp(A,B,Am,Bm,A0,'s');
    % Consegui obter R, S e T corretas
    u(i) = -R(2)*u(i-1)+T(1)*uc(i)-S(1)*y(i)-S(2)*y(i-1);
    % Acho que � esse o bizu
end 

j = 1;
for i=1:N
    if(i<=(j/dt))
        un(i) = u(j);
        ucn(i) = uc(j);
        yn(i) = y(j);
        a1n(i) = a1(j); a2n(i) = a2(j);
        b0n(i) = b0(j); b1n(i) = b1(j);
    else
        j=j+1;
        un(i) = u(j);
        ucn(i) = uc(j);
        yn(i) = y(j);
        a1n(i) = a1(j); a2n(i) = a2(j);
        b0n(i) = b0(j); b1n(i) = b1(j);
    end
end

taux=1:TMAX;

figure (1)
plot(taux,uc,taux,y),legend('u_c','y'),grid;
% Dar zoom no ponto inicial, observa-se o mesmo pico do livro
figure(2)
plot(t,a1n,t,a2n,t,b0n,t,b1n),legend('a_1','a_2','b_0','b_1'),grid;

figure(3)
plot(t,un),grid;

