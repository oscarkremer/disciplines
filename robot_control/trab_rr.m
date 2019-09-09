clear 
clc
close all

m = 1;
t = 0:0.002:24;
t = t';
t_plot = 0:0.002:23.996;
t_plot = t_plot';
sizes = size(t);
%Descomentar para segunda trajetoria
yd = 0.25*sin(t) + 0.25*ones(sizes(1),1);
xd = 0.25*cos(t) + 0.25*ones(sizes(1),1);

q = zeros(sizes(1),2);
c2 = (xd.^2 + yd.^2 - 0.0625*ones(sizes(1),1) - 0.0625*ones(sizes(1),1))/(2*0.25*0.25);
s2 = (1-c2.^2).^0.5;
theta2d = atan2(s2,c2);
k1 = ones(sizes(1),1)*0.25 + 0.25*c2;
k2 = 0.25*s2;
theta1d = atan2(yd,xd) - atan2(k2,k1);

x1d = 0.25*cos(theta1d);

y1d = 0.25*sin(theta1d);


x2d = x1d + 0.25*cos(theta1d + theta2d);
y2d = y1d + 0.25*sin(theta1d + theta2d);

plot(x2d, y2d);