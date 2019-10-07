%%%Clear console, variables and close all windows%%%
clc;
clear ;
close all;
%%%Load data and separate classes%%%
data= load('formantdata.mat');
X = data.D;
Y = data.L;
total_size = size(X); total_size = total_size(1);
folds_number = 7;
size_folds = round(total_size/folds_number);

for i=1:folds_number
  if i==1
    X_test = X(1:size_folds, :);
    X_train = X(size_folds+1:total_size, :);
    Y_test = Y(1:size_folds, :);
    Y_train = Y(size_folds+1:total_size, :);
  elseif i==folds_number
    X_test = X((folds_number-1)*size_folds+1:total_size, :);
    X_train = X(1:(folds_number-1)*size_folds, :);
    Y_test = Y((folds_number-1)*size_folds+1:total_size, :);
    Y_train = Y(1:(folds_number-1)*size_folds, :);
  else
    X_test = X((i-1)*size_folds+1:i*size_folds, :);
    X_train = [X(1:(i-1)*size_folds, :); X(i*size_folds+1:total_size, :)];
    Y_test = Y((i-1)*size_folds+1:i*size_folds, :);
    Y_train = [Y(1:(i-1)*size_folds, :); Y(i*size_folds+1:total_size, :)];
  end
end

x1 = [];
x2 = [];
x1_test = [];
x2_test = [];
sizes_train = size(X_train); sizes_train = sizes_train(1);
sizes_test = size(X_test); sizes_test = sizes_test(1);

for i=1:sizes_train
  if Y_train(i)==0
    x1 = [x1; X_train(i,:)];
  else
    x2 = [x2; X_train(i,:)];
  end
end

for i=1:sizes_test
  if Y_test(i)==0
    x1 = [x1; X_train(i,:)];
  else
    x2 = [x2; X_train(i,:)];
  end
end
%%%Calculate size of each class, mean and covariance matrix%%%

size1 = size(x1);
size1 = size1(1);
mu1 = mean(x1);

size2 = size(x2);
size2 = size2(1);
mu2 = mean(x2);

cov1 = zeros(2,2);
cov2 = zeros(2,2);
for i=1:size1
  cov1(1,1)= cov1(1,1) + (x1(i,1) - mu1(1))^2;
  cov1(1,2)= cov1(1,2) + (x1(i,1) - mu1(1))*(x1(i,2) - mu1(2));
  cov1(2,1)= cov1(2,1) + (x1(i,2) - mu1(2))*(x1(i,1) - mu1(1));
  cov1(2,2)= cov1(2,2) + (x1(i,2) - mu1(2))^2;
end

for i=1:size2
  cov2(1,1)+= (x2(i,1) - mu2(1))^2;
  cov2(1,2)+= (x2(i,1) - mu2(1))*(x2(i,2) - mu2(2));
  cov2(2,1)+= (x2(i,2) - mu2(2))*(x2(i,1) - mu2(1));
  cov2(2,2)+= (x2(i,2) - mu2(2))^2;
end

cov1 = (1/size1)*cov1;
cov2 = (1/size2)*cov2;
cov_t = cov1+cov2;

%%% 3-D Plot of distributions and its contours %%%

x = mvg(mu1', mu2', cov_t);
phi = size2/(size1+size2);
predict_1 = [];
predict_2 = [];

%%% Predict each class from the gaussian distributions generated %%%

for i=1:449
  if pdf(x0(i,:), mu1, cov_t) > pdf(x0(i,:), mu2,cov_t);
     predict_1 = [predict_1; x0(i,:)];
  else
     predict_2 = [predict_2; x0(i,:)];
  end
end

%%% Plots of predictions %%%


figure(3)
clf;
scatter (predict_1(:,1), predict_1(:,2), "r");
hold on;
scatter (predict_2(:,1), predict_2(:,2), "b");
hold on;
mu0 = 0.5*(mu1+mu2);
sigma1 = std(x0(:,1), 1);
sigma2 = std(x0(:,2), 1);
orientation = [(1/sigma2^2)*(mu1(2)-mu2(2)); (1/sigma1^2)*(mu1(1)-mu2(1))];
t = linspace(-2,3,100);
bound = mu0(2) - orientation(2)*(t-mu0(1))/orientation(1);
plot(t,bound);
hold off;
legend ("Classe 1"," Classe 2 ", "Boundary");
grid on;
title ("Predictions - Plot");

%%% Plots of Right Results %%%

figure(4)
clf;
scatter (x1(:,1), x1(:,2), "g");
hold on;
scatter (x2(:,1), x2(:,2), "black");
hold off;
grid on;
legend ("Classe 1"," Classe 2 ");
title ("Original Dataset Classification- Plot");

