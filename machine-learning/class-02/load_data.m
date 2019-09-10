function data = load_datas()
    X = load('data/x.dat');
    y = load('data/y.dat');
    data = [X,y];
end
