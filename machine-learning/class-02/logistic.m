X = load('x.dat');
Y = load('y.dat');
x = zeros(2,1);
res = 50;
for i=1:res,
  for j=1:res,
    x(1) = 2*(i-1)/(res-1) - 1;
    x(2) = 2*(j-1)/(res-1) - 1;
    lambda=0.0001;
    tau=0.5;
    norma =  ones(69,1)*x'-X;
    norma_n = ones(69,1);
    for l=1:69
      norma_n(i) = norm(norma(i,:));      
    end  
    w = exp((-(norma_n).^2)/(2*tau*tau));
    theta  = [0; 0];
    for k=1:100
      grad = X'*w*(Y-((1+exp(-theta'*X')).^-1));
      D = zeros(size(X)(1));
      cal_D = -w*(1/(1+exp(-theta'*X')))*(1-(1/(1+exp(-theta'*X'))));
      for i=1:69
        D(i)(i) = -cal_D(i);
      end
    H = X'*D*X - lambda*eye(2);
    theta=theta-inv(H)*grad;
    end
    y = 1/(1+exp(-theta'*x));
    pred(j,i) = lwlr(X, y, x, tau);
  end
end
figure(1);
clf;
axis off;
hold on;
imagesc(pred, [-0.4 1.3]);
plot((res/2)*(1+X(y==0,1))+0.5, (res/2)*(1+X(y==0,2))+0.5, 'ko');
plot((res/2)*(1+X(y==1,1))+0.5, (res/2)*(1+X(y==1,2))+0.5, 'kx');
axis equal;
axis square;
text(res/2 - res/7, res + res/20, ['tau = ' num2str(tau)], 'FontSize', 18);