pkg load control;
t = linspace(0,10,1000);
x = [];
for i=1:25
  if rand(1,1) > 0.5
    x = [x 5*ones(40,1)'];
  else
    x = [x 5*zeros(40,1)'];
  end
end
plot(t,x);
s = tf('s');
R = 39000;
C = 10^-5;
G = 1/(R*C*s+1);
lsim(G,x,t);
