%%%Clear console, variables and close all windows%%%
clc;
clear ;
close all;
%%%Load data%%%

data= load('formantdata.mat');
x0 = data.D;
y0 = data.L;

%%%Calculate size , mean and covariance matrix%%%

size1 = size(x0);
size1 = size1(1);
mu1 = mean(x0);

x0 = x0 - mu1;

cov1 = zeros(2,2);

for i=1:size1
  cov1(1,1)= cov1(1,1) + (x0(i,1) - mu1(1))^2;
  cov1(1,2)= cov1(1,2) + (x0(i,1) - mu1(1))*(x0(i,2) - mu1(2));
  cov1(2,1)= cov1(2,1) + (x0(i,2) - mu1(2))*(x0(i,1) - mu1(1));
  cov1(2,2)= cov1(2,2) + (x0(i,2) - mu1(2))^2;
end


cov1 = (1/size1)*cov1;

for i=1:size1
  x0(i,:) = (inv((cov1).^0.5)*(x0(i,:)'))';
end

%%%Find eig values, eig vectors and y%%%

[V, lambda] = eig(cov1);
figure(1)
clf;
scatter (x0(:,1), x0(:,2), "r");
hold on;
quiver(0, 0, V(1,1), V(2,1))
quiver(0, 0, V(1,2), V(2,2));
hold off;
grid on;
title ("Principal Component - Plot");


lambda = max(lambda);
if lambda(1)>lambda(2)
  u = V(:,1);
else
  u = V(:,2);
end

y = zeros(size1,1);
for i=1:size1
  y(i) = u'*x0(i,:)';
end


