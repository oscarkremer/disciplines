xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t = 0.:0.01:5
t = t'
S = 3*sin(t) + 2*cos(t) + 6*3*t^2;

PHI = [sin(t) cos(t) 3*t^2] ;
theta1 = inv(PHI(:,1)'*PHI(:,1))*PHI(:,1)'*S;
V1 = 0.5*sum((S - (PHI(:,1)*theta1))^2)

theta2 = inv(PHI(:,1:2)'*PHI(:,1:2))*PHI(:,1:2)'*S;

V2 = 0.5*sum((S - (PHI(:,1:2)*theta2))^2)


theta3 = inv(PHI'*PHI)*PHI'*S;

V3 = 0.5*sum((S - (PHI*theta3))^2)

plot2d(t, PHI(:,1)*theta1,3)
plot2d(t, PHI(:,1:2)*theta2, 4)
plot2d(t, PHI*theta3, 5)yh
plot2d(t, S, 6)

/////////////////////////////////////////////


