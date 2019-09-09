n = linspace(0, 9, 10);

subplot(4,1,1);
k = 1;
y = sin(2*k*pi*n/5);
stem(n,y)
title('Sinal com k=1');
grid on;

subplot(4,1,2); 
k = 2;
y = sin(2*k*pi*n/5);
stem(n,y)
title('Sinal com k=2');
grid on;

subplot(4,1,3); 
k = 4;	
y = sin(2*k*pi*n/5);
stem(n,y)
title('Sinal com k=4');
grid on;

subplot(4,1,4); 
k = 6;
y = sin(2*k*pi*n/5);
stem(n,y);
title('Sinal com k=6');
grid on;