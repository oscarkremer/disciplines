function points = boundary(mean1, mean2, cov1, cov2)
  [x y]=meshgrid(-4:0.1:4,-4:0.1:4);
  x1=[x(:) y(:)]';
  mn=repmat(mean1,1,size(x1,2));
  mulGau1= 1/(2*pi*det(cov1)^(1/2))*exp(-0.5.*(x1-mn)'*inv(cov1)*(x1-mn));
  mn=repmat(mean2,1,size(x1,2));
  mulGau2= 1/(2*pi*det(cov2)^(1/2))*exp(-0.5.*(x1-mn)'*inv(cov2)*(x1-mn));
  G1=reshape(diag(mulGau1),81,81);
  G2=reshape(diag(mulGau2),81,81);
  x_out = [];
  y_out = [];
  for i=1:81
    for j=1:81
      if abs(G1(i,j)- G2(i,j)) < 0.001 && abs(G1) > 0.00001 && abs(G2) > 0.00001

        x_out = [x_out; x(i,j)];
        y_out = [y_out; y(i,j)];
      end
    end 
  end 
  points = [x_out, y_out];
 end

