pkg load control

A = [1 0;0 -2];
B = [1 1]';
C = [1 1];
D = 0;
L = 1/3*[16 -1]'
K = [-2 0];
I = eye(2);
Aa = [A B*K; L*C (A-L*C+B*K)];
Ba = [B; B];
Ca = [C zeros(size(C))];
t = 0:0.01:5;
u = 0*t;
x0 = [1 1 0 0]'
[Y,X] = lsim(Aa, Ba, Ca, D, u,t, x0);
E = [X(:,1)-X(:,3) X(:,2)-X(:,4)];

