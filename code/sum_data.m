%СОПРОТИВЛЕНИЕ АНТЕННЫ
% ANTRES = sum(sum(ZAR)); 

%%Таблица сопротивлений (I;J) вибратора  в виде
% N (1..m*n), dx/lambda, dy/lambda, R, X
% N (m*n..2*m*n), dx/lambda, dy/lambda, R, X  (Для вибраторов рефлектора)
I = 3; J=2;  %координаты вибратора (строка;столбец)
%%


Excl(m*n, 5) = 0;
Excl(1:m*n,1) = 1:18;
Excl(m*n+1:2*m*n,1) = 1:18;

s=0; sum1=0; sum2=0;
for i = 1:6
    for j = 1:3
    s=s+1;
    Excl(s,2) = DH(I,J,s);
    Excl(s+18,2) = DH(I,J,s+18);
    Excl(s,3) = DY(I,J,s);
    Excl(s+18,3) = DY(I,J,s+18);
    Excl(s,4) = real(ZA(I,J,s));
    Excl(s+18,4) = real(ZA(I,J,s+18));
    sum1  = sum1 + ZA(I,J,s);
    sum2  = sum2 + ZA(I,J,s+18);
    Excl(s,5) = imag(ZA(I,J,s));
    Excl(s+18,5) = imag(ZA(I,J,s+18));
     
    end
end
%% Полное сопротивление (I;J) вибратора 
VibRes = sum1 + sum2*(-1);



