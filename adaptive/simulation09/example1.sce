xdel(winsid());//fechar janelas de figuras
clear();//limpar memória
clc();//limpar console
t0 = 0;
delta_t = 0.001;
final_time = 10;
t = (0:delta_t:final_time);
t = t'
ifinal=size(t);ifinal=ifinal(1);
sampling_rate = 0.5
ratio = sampling_rate/delta_t;
y = zeros(ifinal);
u = zeros(ifinal);



entradac = 2*ones(ifinal,1);
for i=1:ifinal
     if (i-1 - ratio*fix((i-1)/ratio)==0) then
        if (i==1) then
            y_2 = 0;
            y_1 = 0;
            u_1 = 0;
         elseif i==2 then
            y_1 = y;
            y_2 = 0;
            u = entradac(i-1);       
     end  
       
             
         
end



