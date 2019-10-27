%%%% faz o plot tridimensional das distribuicoes e dos contornos das gaussianas %%%%

function x=mvg(mean1, mean2, cov1, cov2)
  [x y]=meshgrid(-4:0.1:4,-4:0.1:4);
  x1=[x(:) y(:)]';
  mn=repmat(mean1,1,size(x1,2));
  mulGau1= 1/(2*pi*det(cov1)^(1/2))*exp(-0.5.*(x1-mn)'*inv(cov2)*(x1-mn));
  mn=repmat(mean2,1,size(x1,2));
  mulGau2= 1/(2*pi*det(cov2)^(1/2))*exp(-0.5.*(x1-mn)'*inv(cov2)*(x1-mn));
  G1=reshape(diag(mulGau1),81,81);
  G2=reshape(diag(mulGau2),81,81);
  figure(1);
  clf;
  mesh(-4:0.1:4, -4:0.1:4, G1);
  hold on;
  mesh(-4:0.1:4, -4:0.1:4, G2);
  legend('Classe 1 - y=0', 'Classe 2 -y=1');
  hold off;
  figure(2);
  clf;
  contour(-4:0.1:4, -4:0.1:4, G1);
  hold on;
  contour(-4:0.1:4, -4:0.1:4, G2);
  hold off;
  grid on;
end

