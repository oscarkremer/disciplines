%question 02 matlab list

a = [1];
b = [1; 1];
wmax = 1000;
w = linspace(0, 2*pi, 1000)';

H = freqresp(b, a, w);

figure(1)
plot (w, 20*log10(abs(H)))

xlim ([0 2*pi])
ylim ([-10 10])

figure(2)
plot (w, 180*angle(H)/pi)

xlim ([0 2*pi])
ylim ([-200 200])


a = [1];
b = [1; -1];
H = freqresp(b, a, w);
figure(3)
plot (w, 20*log10(abs(H)))
xlim ([0 2*pi])
ylim ([-10 10])
figure(4)
plot (w, 180*angle(H)/pi)
xlim ([0 2*pi])
ylim ([-200 200])
