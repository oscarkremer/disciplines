
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
mi1 = 0;
mi2 = 0;
mi = [mi1;mi2];
sigma = [1 0; 0.2 5];

x = mvg(mi, sigma);