clear
clc
close all

delta = [1; zeros(99,1)];
a = [1; 3/5; zeros(98,1)];
b = [1; zeros(99,1)];
h1 = zeros(20,1);
h2 = zeros(20,1);

for i=0:40
  if i==0
    h1(i+1) = 1;
    h2(i+1) = 1;
  else
    h1(i+1) = 0.6*h1(i);
    h2(i+1) = (0.6^i)*h2(i);
  end
end

figure(1)
stem(linspace(0,40,41), h1)
title('Primeiro Sinal - h1')


figure(2)
stem(linspace(0,40,41), h2)
title('Segundo Sinal - h1')


z1 = conv(h1, ones(41,1));
z2 = conv(h2, ones(41,1));

s1 = zeros(79,1);
s2 = zeros(79,1);

for i=0:78
  if i==0
    s1(i+1) = 1;
    s2(i+1) = 1;
  else
    s1(i+1) = 0.6*s1(i) + 1;
    s2(i+1) = (0.6^i)*s2(i) + 1;
  end
end

figure(3)
subplot(2,1,1);
stem(linspace(0, 80, 81), z1);
title('Primeiro Sinal - z1')

subplot(2,1,2); 
stem(linspace(0, 78, 79), s1);
title('Primeiro Sinal - s1')

figure(4)
title('Segundo Sinal')
subplot(2,1,1);
stem(linspace(0, 80, 81), z2);
title('Segundo Sinal - z2')

subplot(2,1,2); 
stem(linspace(0, 78, 79), s2);
title('Segundo Sinal - s2')
