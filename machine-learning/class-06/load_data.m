function data = load_data()
    X = load('x.dat');
    y = load('y.dat');
    data = [X,y];
end
