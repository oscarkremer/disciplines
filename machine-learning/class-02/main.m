fdata = load_data();
X = data(:,1:2);
y = data(:,3);

plot_lwlr(X, y, 0.01, 40);
plot_lwlr(X, y, 0.05, 40);
plot_lwlr(X, y, 0.1, 40);
plot_lwlr(X, y, 0.5, 40);
plot_lwlr(X, y, 1, 40);
plot_lwlr(X, y, 5, 40);
