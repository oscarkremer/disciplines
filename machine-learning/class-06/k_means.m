%%%Clear console, variables and close all windows%%%
clc;
clear ;
close all;

%%% Carregar dados e inicializacao das variaveis %%%
data= load('formantdata.mat');
x0 = data.D;
y0 = data.L;
x1 = [];
x2 = [];
mu1 = [randn randn];
mu2 = [randn randn];
epsilon = 0.001;
max = 10;
J_old = 100000;
J = zeros(max,1);

%%% Definir o laco principal que pode ser executado um numero maximo de vezes %%%

for j=1:max
  m = size(x0); m = m(1);
  class1 = [];
  class2 = [];
  %%% Calcular as distancia com base nos valores dos centros e efetuar a classificacao %%% 
  for i=1:m
    dist1 = (x0(i,1) - mu1(1))^2 + (x0(i,2) - mu1(2))^2;
    dist2 = (x0(i,1) - mu2(1))^2 + (x0(i,2) - mu2(2))^2;
    if dist1 < dist2
      class1 = [class1; x0(i,:)];
    else
      class2 = [class2; x0(i,:)];
    end
  end
  
  %%% Plot de dispersao dos dados para o passo j %%% 

  figure(j)
  clf;
  scatter (class1(:,1), class1(:,2), "g");
  hold on;
  scatter (class2(:,1), class2(:,2), "black");
  %%% Demarcacao dos centros de cada cluster com -x %%% 
  plot(mu1(1),mu1(2),'-x');
  plot(mu2(1),mu2(2),'-x');
  hold off;
  grid on;
  legend ("Classe 1"," Classe 2 ");
  title ("Original Dataset Classification- Plot");
  %%% Calcular a distorcao total dos dados %%% 
  J_new = distortion(mu1, mu2, class1, class2);
  mu1 = mean(class1);
  mu2 = mean(class2);
  J(j) = J_new;
  %%% Verificar se o criterio de parada e atendido %%% 
  if (1-J_new/J_old) < epsilon
    break;
  else
    J_old = J_new;
  end
end

%%% Plot da distorcao ao longo dos passos de treinamento %%% 

figure(11)
clf;
plot(linspace(1,10,10), J);
axis([1 max  0 1000]);
grid on;
title ("Distortion through Steps");