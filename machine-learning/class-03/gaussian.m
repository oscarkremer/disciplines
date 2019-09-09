clear 

function x=mvg(mean,cov)
  %%%this is the xy plane limits
  [x y]=meshgrid(-4:0.1:4,-4:0.1:4);
  x1=[x(:) y(:)]';
  %%%multivar Gassiaan
  mn=repmat(mean,1,size(x1,2));
  mulGau= 1/(2*pi*det(cov)^(1/2))*exp(-0.5.*(x1-mn)'*inv(cov)*(x1-mn));
  G=reshape(diag(mulGau),81,81);
  figure(1);
  mesh(-4:0.1:4, -4:0.1:4, G);
  figure(2);
  contour(-4:0.1:4, -4:0.1:4, G);
  %plot3(x1(1,:),x1(2,:),diag(mulGau))
endfunction

data= load('formantdata.mat');
x = data.D;
y = data.L;
x0 = [];
x1 = [];
for i=1:449
  if y(i)==0
    x0 = [x0; x(i,:)];
  else
    x1 = [x1; x(i,:)];
  end
end

mu0 = mean(x0);
mu1 = mean(x1);
cov = zeros(2,2);
sum11 = 0;
sum12 = 0;
sum21 = 0;
sum22 = 0;

for i=1:262
  sum11 = sum11 + (x0(i,:) - mu0)'*(x0(i,:) - mu0);
end
for i=1:187
  sum12 = sum22 + (x1(i,:) - mu1)'*(x0(i,:) - mu0);
  sum21 = sum22 + (x0(i,:) - mu0)'*(x1(i,:) - mu1);
  sum22 = sum22 + (x1(i,:) - mu1)'*(x1(i,:) - mu1);
end
cov = sum11+sum12+sum21+sum22;
cov = (1/449)*cov;
x = mvg(mu0', cov);
x = mvg(mu1', cov);
phi = 187/449;
