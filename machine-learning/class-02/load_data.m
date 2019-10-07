function data = load_data()
    X = load('data/x.dat');
    y = load('data/y.dat');
    data = [X,y];
end
