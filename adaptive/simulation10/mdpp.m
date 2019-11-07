% ROTINA PARA IMPLEMENTAÇÃO DO ALGORITMO 3.1
function [R,S,T] = mdpp(A,B,Am,Bm,A0,op)

% Condições de compatibilidade
%length(A) = length(Am)
%length(B) = length(Bm)
%length(A0) = length(A)-length(Bplus)-1
%length(A0) = d0-1
%Bm = conv(Bminus,Bml)

if op == 'c'
    % Passo 1
    Bplus = B/B(1); % Tornando B um polinômio mônico
    [Bminus,rBminus] = deconv(B,Bplus);
    [Bml, rBml] = deconv(Bm,Bminus);

    % Passo 2
	A0Am = conv(A0,Am);
    [Rl,BminusS] = deconv(A0Am,A); 
    
    R = conv(Rl,Bplus);
    S = deconv(BminusS,Bminus);
    if S(1)==0
        S = S(2:end);
    end
    T = conv(A0,Bml);
end

if op == 's'
    % Passo 1
    Bplus = 1; 
    Bminus = B;
    beta = polyval(Am,1)/polyval(B,1);
    %Bm = beta*B;
    T = beta*A0;
    
    % Passo 2
	Ac = conv(A0,Am);
    
    % B tem que ser do msm tamanho de A
    aux = B;
    B = zeros(1,length(A));
    B(2:end) = aux;
    
    t = Ac';

    Sab(:,1) = [A 0];
    Sab(:,2) = [0 A];
    Sab(:,3) = [B 0];
    Sab(:,4) = [0 B];
    
    w = inv(Sab)*t;
    w = w';
    
    R = w(1:2);
    S = w(3:4);
    %[R,BminusS] = deconv(Ac,A); 
    
    %R = conv(Rlinha,Bplus);
    %S = deconv(BminusS,Bminus);
    %if S(1)==0
    %    S = S(2:end);
    %end
    
end

end

