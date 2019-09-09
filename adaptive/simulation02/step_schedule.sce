xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console

function gain=g(amp)
    amp1 = 0.5;
    amp2 = 1.5;
    amp3 = 2.5;
    kl2 = 4*amp2^(3/4);
    kl3 = 4*amp3^(3/4);
    if amp <= 0.5 then
        k = 0.10;
        kl = 4*amp1^(3/4);
    elseif amp <= 1.5 then
        k = 0.025;
        kl = 4*amp2^(3/4);    
    else
        k = 0.012575;
        kl = 4*amp3^(3/4);        
    end
    gain = [k, kl];
endfunction


gain = g(0.3)
gain(1)

num=poly([gain(1)*gain(2)],'s','coeff');
den=poly([gain(1)*gain(2) 1 2 1],'s','coeff');
g1=syslin('c',num/den) //define tf
gs1=csim(0.3*ones(1, 1001), 0:0.05:50, g1);

gain = g(1.0)
num = poly([gain(1)*gain(2)],'s','coeff');
den = poly([gain(1)*gain(2) 1 2 1],'s','coeff');
g2 = syslin('c',num/den) //define tf
gs2=csim(1.0*ones(1, 1001), 0:0.05:50, g2);

gain = g(2.0)
num = poly([gain(1)*gain(2)],'s','coeff');
den = poly([gain(1)*gain(2) 1 2 1],'s','coeff');
g3 = syslin('c',num/den) //define tf
gs3=csim(2.0*ones(1, 1001), 0:0.05:50, g3);

figure(0)
plot2d(0:0.05:50, gs1); // plotting step response
plot2d(0:0.05:50, gs2); // plotting step response
plot2d(0:0.05:50, gs3); // plotting step response

figure(1) 
gs2_complete = gs2+0.3*ones(1,1001);
gs3_complete = gs3+1.3*ones(1,1001);

plot2d([0:0.05:50 50:0.05:100 100:0.05:150] , [gs1 gs2_complete gs3_complete]); // plotting step response



