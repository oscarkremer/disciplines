
clear; close all; clc; 
%% Definicao dos valores de treino e normalização dos mesmos
rooms = randi([1 5],1,500)/5;
area = randi([100 1000],1, 500)/1000;
rooms = sort(rooms);
area = sort(area);
price = randi([50, 100],1, 500)/100;
price = sort(price);
%% Definição do learning rate e dos valores iniciais para os coeficientes
alpha= 0.001;
theta0i = 0;
theta1i = 0;
theta2i = 0;
theta0 = zeros(1,100);
theta1 = zeros(1,100);
theta2 = zeros(1,100);
n = 100;
%% Laço para treinamento
for k = 1:n
    h = theta0i*ones(1,500) + theta1i*area + theta2i*rooms;
    theta0f = theta0i + alpha*(price - h)*ones(1,500)';
    theta0i = theta0f;
    h = theta0i*ones(1,500) + theta1i*area + theta2i*rooms;
    theta1f = theta1i + alpha*(price - h)*rooms';
    theta1i = theta1f;
    h = theta0i*ones(1,500) + theta1i*area + theta2i*rooms;
    theta2f = theta2i + alpha*(price - h)*area';
    theta2i = theta2f; 
    theta0(k) = theta0f;
    theta1(k) = theta1f;
    theta2(k) = theta2f;
end

%% Laço para plot dos pesos ao longo do período de treinamento
figure(1)
plot(theta0)
title('Theta 0')
figure(2)
plot(theta1)
title('Theta 1')
figure(3)
plot(theta2)
title('Theta 2')
%%% Plot da dispersão dos dados e plot da reta gerada pelo LMS para o caso
%%% específico
figure(4)
x_axis = linspace(100,1000,500)/1000;
rooms_axis  = linspace(1,5,500)/5;
h = theta0i*ones(1,500) + theta1i*x_axis + theta2i*rooms_axis;
scatter(1000*area, 100*price,'filled')
hold on;
plot(1000*x_axis, 100*h);

