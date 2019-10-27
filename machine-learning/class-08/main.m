close all; clear;
d = 2;
%processamento dos dados
data= load('formantdata.mat');
X = data.D';
label = data.L';
n = size(X); n=n(2)


% classificacao dos dados

m = floor(n/2);
[z1,model,llh] = mixGaussEm(X);

plotClass(X,z1,1);

x1 = [];
x2 = [];
for i=1:n
  if z1(i)==1
    x1 = [x1; X(:,i)'];
  else
    x2 = [x2; X(:,i)'];
  end
end

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

%%% Plot 3-D e contornos das gaussianas %%%

x = mvg(mu1', mu2', cov1, cov2);
