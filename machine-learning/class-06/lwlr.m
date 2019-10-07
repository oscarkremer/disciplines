function y = lwlr(X_train, y_train, x, tau)
    size_train = size(X_train);
    size_train = size_train(1);
    weights = zeros(size_train,1);
    
    lambda = 0.0001;
    for i=1:size_train
        weights(i) = exp(-((norm(x-X_train(i,:)))^2)/(2*tau*tau));
    end
    theta = [20; 10];
    for j=1:50
        h_theta = 1./(1+exp(-theta'*X_train'))';
        D = eye(size_train);
        aux =   -weights.*h_theta.*(ones(size_train,1)-h_theta);
        for i=1:size_train
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
    
    