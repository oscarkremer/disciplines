%%%Clear console, variables and close all windows%%%
clc;
clear all;
close all;
%%%Load data and separate classes%%%
data= load('formantdata.mat');
X = data.D;
Y = data.L;
%% calculo da media dos dados
mu_total = mean(X);

%% separar as duas classes presentes
x1 = X(find(Y==0),:);
x2 = X(find(Y==1),:);

%% chamar funcao para encontraro autovetor
[D, W_lda] = lda(X, Y); 

Xm = bsxfun(@minus, X, mean(X));

%% fazer a representacao dos dados reduzidos no plano XY
z = Xm*W_lda(:,1);
% and reconstruct it
p = z*W_lda(:,1)';
p = bsxfun(@plus, p, mu_total);
y1 = p(find(Y==0),:);
z1 = z(find(Y==0),:);

y2 = p(find(Y==1),:);
z2 = z(find(Y==1),:);
scale = 5;
figure(1)
pc1 = line([mu_total(1) - scale * W_lda(1,1) mu_total(1) + scale*W_lda(1,1)], [mu_total(2) - scale*W_lda(2,1) mu_total(2) + scale*W_lda(2,1)]);
pc1 = line([mu_total(1) - scale * W_lda(1,2) mu_total(1) + scale*W_lda(1,2)], [mu_total(2) - scale*W_lda(2,2) mu_total(2) + scale*W_lda(2,2)]);
hold on;
set(pc1, 'color', [1 0 0], "linestyle", "--");
p1 = plot(y1(:,1),y1(:,2),"ro", "markersize", 10, "linewidth", 1); 
p2 = plot(y2(:,1), y2(:,2),"go", "markersize", 10, "linewidth", 1);
grid


%% Criar gaussianas para cada uma das distribuicoes

axis1  = 1.25*min(z1):0.001:1.25*max(z1);
pdf1 = mvnpdf(axis1', mean(z1), std(z1));
axis2  = 1.25*min(z2):0.001:1.25*max(z2);
pdf2 = mvnpdf(axis2', mean(z2), std(z2));

%% Fazer plot das gaussianas junto com os valores no eixo

figure(2)
plot(axis1,pdf1,'r');hold on;
plot(axis2,pdf2, 'k');
plot(z1, zeros(size(z1)),'rx');
plot(z2, zeros(size(z2)),'kx');
grid ;

legend ("Classe 1"," Classe 2 ");
title ("Original Dataset Classification- Plot");