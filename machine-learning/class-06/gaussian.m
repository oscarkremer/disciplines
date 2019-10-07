%%%Clear console, variables and close all windows%%%
clc;
clear ;
close all;
%%%Load data and separate classes %%%
data= load('formantdata.mat');
X = data.D;
Y = data.L;
total_size = size(X); total_size = total_size(1);
folds_number = 7;
size_folds = round(total_size/folds_number);
error = zeros(folds_number, 1);
for k=1:folds_number
  if k==1
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
    X_test = X((k-1)*size_folds+1:k*size_folds, :);
    X_train = [X(1:(k-1)*size_folds, :); X(k*size_folds+1:total_size, :)];
    Y_test = Y((k-1)*size_folds+1:k*size_folds, :);
    Y_train = [Y(1:(k-1)*size_folds, :); Y(k*size_folds+1:total_size, :)];
  end
  x1 = [];
  x2 = [];
  x1_test = [];
  x2_test = [];
  size_train = size(X_train); size_train = size_train(1);
  size_test = size(X_test); size_test = size_test(1);

  for i=1:size_train
    if Y_train(i)==0
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
  right = 0;
  wrong = 0;
  for i=1:size_test
    if pdf(X_test(i,:), mu1, cov_t) > pdf(X_test(i,:), mu2,cov_t);
       if Y_test(i) == 0
        right += 1;
       else
        wrong += 1;
       end
    else
       if Y_test(i) == 1
        right += 1;
       else
        wrong += 1;
       end
    end
  end
  error(k) = wrong/(right+wrong);
end

error_medio_percentual = mean(error)*100 