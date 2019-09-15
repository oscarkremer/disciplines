%%% Codigo feito para a disciplina de controle multivariável
%%% Oscar Schmitt Kremer - 2019 - Segundo Semestre
%%% 14 de setembro de 2019
%%% Professor Doutor Carlos Mendes Richter
%% Limpar console, variáveis e fechar janelas
clear;
close all;
clc;
%% definir numero de pontos no gráfico e valores máximo e mínimo de frequência utilizados

number_points = 50000;
down_omega= -100;
upp_omega = 100;

%% Definir o dominio s e a matriz de transferência G

omega = linspace(down_omega, upp_omega,number_points);
s_omega = i*omega;
s = tf('s');
G =  (1/(1.25*(s+1)*(s+2)))*[s-1 s;-6 s-2];

%% Calcular os autovalores

lambda1 = [];
lambda2 = [];
for k=1:number_points
    G_ = evalfr(G,s_omega(k));
    values = eig(G_);
    lambda1 = [lambda1; values(1)];
    lambda2 = [lambda2; values(2)];
end


%% Fazer plot dos characteristic loc separados e juntos

figure(1);
plot(lambda1,'-');
title('Characteristic loci 1')
grid on;

figure(2);
plot(lambda2,'-');
title('Characteristic loci 2')
grid on;

figure(3);
hold on;
plot(lambda2,'-');
plot(lambda1,'-');
pbaspect([1 1 1]);
title('Nyquist for MIMO System')
grid on;