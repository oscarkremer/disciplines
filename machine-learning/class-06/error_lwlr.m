function error = error_lwlr(X_train, y_train, X_test, Y_test, tau, res)
    size_test = size(X_test);
    size_test = size_test(1);        
    pred = zeros(size_test,1);
    right = 0;
    wrong = 0;
    for i=1:size_test
      pred(i) = lwlr(X_train, y_train, X_test(i,:), tau); 
      if pred(i)==Y_test(i)
        right +=1;
      else
        wrong+=1;
      end
    end
    error = wrong/(right+wrong);
end
