data = load_data();
X = data(:,1:2);
Y = data(:,3);
total_size = size(X); total_size = total_size(1);
folds_number = 7;
size_folds = round(total_size/folds_number);
error = zeros(folds_number, 1);
tau = 0.1;
for i=1:folds_number
  if i==1
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
    X_test = X((i-1)*size_folds+1:i*size_folds, :);
    X_train = [X(1:(i-1)*size_folds, :); X(i*size_folds+1:total_size, :)];
    Y_test = Y((i-1)*size_folds+1:i*size_folds, :);
    Y_train = [Y(1:(i-1)*size_folds, :); Y(i*size_folds+1:total_size, :)];
  end
  error(i) = error_lwlr(X_train, Y_train, X_test, Y_test, tau);
end
mean(error)*100