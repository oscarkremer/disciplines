pkg load control
pkg load signal

num = {[1/1.25, -1/1.25], [1/1.25 0]; [-6/1.25], [1/1.25, -2/1.25]};
den = {[1, 3, 2], [1, 3, 2]; [1, 3, 2], [1, 3, 2]};
[a,b,c,d] = tf2ss(num,den)

