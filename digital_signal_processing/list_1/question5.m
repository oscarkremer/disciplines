clear
clc
close all


function y = diffeqn(a, x, yn1)
  y = zeros(31,1);
  sizes = size(x);
  N = sizes(1);
  for i=1:N
    if i==1
      y(i) = a*yn1 + x(1);
    else
      y(i) = a*y(i-1) + x(i);
    end
  end
endfunction

x1 = [1;zeros(30,1)];
y1 = diffeqn(1,x1,0);
x2 = ones(31,1);
y2 = diffeqn(1,x2,0);

figure(1)
clf;
stem(linspace(0,30,31), y1);
ylim([0,1.1]);
title('Resposta 5.1');

figure(2)
stem(linspace(0,30,31), y2);
title('Resposta 5.2');


x1 = ones(31,1);
y1 = diffeqn(1,x1,-1);
x2 = 2*x1;
y2 = diffeqn(1,x2,-1);



figure(3)

n = linspace(0, 9, 10);

subplot(3,1,1);
stem(linspace(0,30,31), y1);
title('Resposta y1');

subplot(3,1,2); 
stem(linspace(0,30,31), y2);
title('Resposta y2');

subplot(3,1,3); 
stem(linspace(0,30,31), 2*y1-y2);
title('Resposta 2y1-y1');

x1 = ones(31,1);
y1 = diffeqn(1,x1,-1);
x2 = 2*x1;
y2 = diffeqn(1,x2,-1);



figure(4)

x = ones(31,1);
y1 = diffeqn(0.5,x,0);
y2 = diffeqn(0.5,x,3);

subplot(2,1,1);
stem(linspace(0,30,31), y1);
title('Resposta y1');

subplot(2,1,2); 

stem(linspace(0,30,31), y2);
title('Resposta y2');

