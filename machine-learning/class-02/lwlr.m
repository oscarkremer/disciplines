function y = lwlr(X_train, y_train, x, tau)
    sizes = size(X_train);
    weights = zeros(sizes(1),1);
    lambda = 0.0001;
    for i=1:sizes(1)
        weights(i) = exp(-((norm(x-X_train(i,:)))^2)/(2*tau*tau));
    end
    theta = [20; 10];
    for j=1:10
        h_theta = 1./(1+exp(-theta'*X_train'))';
        D = eye(69);
        aux =   -weights.*h_theta.*(ones(69,1)-h_theta);
        for i=1:69
          D(i,i) = aux(i);
        end
        z = weights.*(y_train - h_theta);
        grad_theta = X_train'*z - lambda*theta;
        H = X_train'*D*X_train - lambda*eye(2);
        theta = theta - inv(H)*grad_theta;
    end
    if 1./(1+exp(-theta'*x')) > 0.5
        y = 1;
    else
        y = 0;
    end
end
    
    