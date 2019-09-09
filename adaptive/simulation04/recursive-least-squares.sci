xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t = 0.:0.01:5'
S = 3*sin(t) + 2*cos(t) + 6*3*t.^2;
P = 10000*eye(3, 3);
PHI = [sin(t);cos(t);3*t.^2] ;

theta = [0;0;0];
for i=2:501
    K = P*PHI*inv(1+PHI'*P*PHI);
    P = (eye(3,3) - K*PHI')*P;
    theta = theta + K*(S - PHI'*t);
end


//plot2d(t, PHI(:,1)*theta1,3)
//plot2d(t, PHI(:,1:2)*theta2, 4)
//plot2d(t, PHI*theta3, 5)yh
//plot2d(t, S, 6)

/////////////////////////////////////////////


