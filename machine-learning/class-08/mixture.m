%%%Clear console, variables and close all windows%%%
clc;
clear ;
close all;

%%% Carregar dados e inicializacao das variaveis %%%
data= load('formantdata.mat');
x0 = data.D;
y0 = data.L;
mu1 = [-1 1];
mu2 = [1 -1];
cov1 = [rand 0; 0 rand];
cov2 = [rand 0; 0 rand];

max = 10;
m = size(x0); m = m(1);

weight1 = zeros(m, 1);
weight2 = zeros(m, 1);
x1 = [];
x2 = [];
for i=1:449
  if y0(i)==0
    x1 = [x1; x0(i,:)];
  else
    x2 = [x2; x0(i,:)];
  end
end
%%% Definir o laco principal que pode ser executado um numero maximo de vezes %%%

for j=1:15
  %%% Calcular as distancia com base nos valores dos centros e efetuar a classificacao %%% 
  for i=1:m
    weight1(i) = pdf(x0(i,:), mu1, cov1);
    weight2(i) = pdf(x0(i,:), mu2, cov2);
  end
  cov1 = zeros(2,2);
  cov2 = zeros(2,2);
  mu1(1) = sum(weight1.*x0(:,1))/sum(weight1);
  mu1(2) = sum(weight1.*x0(:,2))/sum(weight1);
  mu2(1) = sum(weight2.*x0(:,1))/sum(weight2);
  mu2(2) = sum(weight2.*x0(:,2))/sum(weight2);

  for i=1:m
    cov1(1,1)= cov1(1,1) + weight1(i)*(x0(i,1) - mu1(1))^2;
    cov1(1,2)= cov1(1,2) + weight1(i)*(x0(i,1) - mu1(1))*(x0(i,2) - mu1(2));
    cov1(2,1)= cov1(2,1) + weight1(i)*(x0(i,2) - mu1(2))*(x0(i,1) - mu1(1));
    cov1(2,2)= cov1(2,2) + weight1(i)*(x0(i,2) - mu1(2))^2;  
    cov2(1,1)= cov2(1,1) + weight2(i)*(x0(i,1) - mu2(1))^2;
    cov2(1,2)= cov2(1,2) + weight2(i)*(x0(i,1) - mu2(1))*(x0(i,2) - mu2(2));
    cov2(2,1)= cov2(2,1) + weight2(i)*(x0(i,2) - mu2(2))*(x0(i,1) - mu2(1));
    cov2(2,2)= cov2(2,2) + weight2(i)*(x0(i,2) - mu2(2))^2;  
  end
  cov1 = (1/sum(weight1))*cov1;
  cov2 = (1/sum(weight2))*cov2;
  
end

x = mvg(mu1', mu2', cov1, cov2);
figure(4)
clf;
scatter (x1(:,1), x1(:,2), "g");
hold on;
scatter (x2(:,1), x2(:,2), "black");
hold off;
grid on;
legend (" Classe 1"," Classe 2 ");
title ("Original Dataset Classification- Plot");

