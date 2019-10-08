close all;
%%

Freq.lambda0 = 0.1968:-0.001:0.1637;        %Предполагаемый диапазон длин работы антенны
Freq.lambda0_= 0.18;        %Рабочая длина волны


Freq.lambda0 = 0.21:-0.001:0.15;   

Freq.C=299792458;  
Freq.F = Freq.C./Freq.lambda0; 

%%
CD.L = 0.09; B.L = 0.06; A.L = 0.07;        %Длины кабеля в секциях
CD.ky = 1.19; B.ky =1.52; A.ky = 1.22;      %Коэффициенты укорочения
CD.W = 150;  B.W=75;  A.W=50;               %волновые сопротивления
CD.Wcy = 65.84; CD.lcy =Freq.lambda0_/4;
%Длина волны в коакс. фидерах
CD.lamCD =Freq.lambda0./CD.ky;
B.lamB =Freq.lambda0./B.ky;
A.lamA =Freq.lambda0./A.ky;

% %Волновые споротивления четвертьволнового трансформатора 
B.Wtr = 37.5; A.Wtr=37.5; 
B.Ltr = (Freq.lambda0_/B.ky)/4;  %Длины четрвертьволновых трансформаторов
A.Ltr = (Freq.lambda0_/B.ky)/4;
% %Волновые споротивления трехступенчатых трансформаторов
% B.W1 = 43.63; B.W2 = 33.23; B.W3 = 25.31;    

% A.W1 = 44.54; A.W2 = 35.35; A.W3 = 28.06;
% Freq.L = 0.1513/4;    % Длина секции трехступенчатого трансформатора

% Freq.beta_coax(k) =6.28./Freq.lam1(k);
%%
CD.beta_coax =6.28./CD.lamCD;
B.beta_coax =6.28./B.lamB;
A.beta_coax =6.28./A.lamA;



 for k=1:length(Freq.lambda0) 
  
%     k = 1;                                            
                                              
%%
                                             lam = Freq.lambda0(k);                                            

%%
% Пересчет сопротивлений решетки
ZAR(m,n) = 0;       %Наведенные сопротивления вибраторов решетки
ZA(m,n,m*n*2) = 0;  % Массив сопротивлений R,jX. ...
                    %(i,j, 1..m*n) - взаимные сопротивления i,j вибраторов антенной решетки
                    %(i,j, m*n..2*m*n) - взаимные сопротивления i,j вибратора и вибраторов
                    % рефлектора
                    
for ii =1:m     % 1:(round(m/2))  
for jj=1:n          %1:(round(n/2))
S=0;
for i0 =1:m
for j0=1:n
S=S+1; 
%% EMF

if (ii~=i0)&&(jj~=j0)
ZA(ii,jj,S)=emf(x,y, DY(ii,jj,S),lam, DH(ii,jj,S));  % Z от других вибраторов антенны  
elseif (ii==i0)&&(jj==j0)  %собственное сопротивление вибратора
ZA(ii,jj,S) = Zvib;
else 
ZA(ii,jj,S)=emf(x,y, DY(ii,jj,S),lam, DH(ii,jj,S));  
end
ZA(ii,jj,S+m*n)=emf(x,y, DY(ii,jj,S+m*n),lam, DH(ii,jj,S+m*n));  %Z от рефлектора
%% 

%%Суммирование наведенных сопротивлений вибратора от других вибраторов
%%и рефлектора
 ZAR(ii,jj) = ZAR(ii,jj)+ ZA(ii,jj,S) +  (ZA(ii,jj,S+m*n)*(-1));

end
end

end

end
clear  ii i0 jj j0;
clear hh bb dd ZA S;

%%     
                                           clear lam;

%%

                              %Блок C и D                                 
for i1 =1:3   
    %Сопротивления с учетом согласующего стакана
    
%     Z.C(i1) = 1/(( 1/(j*CD.Wcy*tan(6.28*CD.lcy/Freq.lambda0(k))) + 1/ZAR(i1,1)));
  Z.C(i1)=ZAR(i1,1);
    Z.C(i1) = CD.W *  (  Z.C(i1)+ j*CD.W*tan(CD.beta_coax(k)*CD.L) )/(CD.W+Z.C(i1)*j*tan(CD.beta_coax(k)*CD.L));
    
%     Z.C(i1) = W* (ZAR(i1,1)+j*W*tan(Freq.beta_coax(k)*CD.L))/(W+j*ZAR(i1,1)*tan(Freq.beta_coax(k)*CD.L));  %Трансформация сопротивления в линии
%     Z.C(i1) =1/ ((1/(Z.C(i1))) + (1/ (j*W*tan(Freq.beta_coax(k)*(Freq.lam1_0/4)))));  %Влияние четвертьволнового стакана
   
%     Z.D(i1) = (W* (ZAR(i1,2)+j*W*tan(Freq.beta_coax(k)*CD.L))/(W+j*ZAR(i1,2)*tan(Freq.beta_coax(k)*CD.L)));
%     Z.D(i1) = 1/ ((1/(Z.D(i1))) + (1/ (j*W*tan(Freq.beta_coax(k)*(Freq.lam1_0/4)))));
%     Z.D(i1) = 1/( 1/(j *CD.Wcy*tan(6.28*CD.lcy/Freq.lambda0(k)))+1/ZAR(i1,2));
  Z.D(i1)=ZAR(i1,2);    
Z.D(i1) = CD.W * (Z.D(i1) + j*CD.W*tan(CD.beta_coax(k)*CD.L))/(CD.W+Z.D(i1)*j*tan(CD.beta_coax(k)*CD.L));
    
    
    Ksv.VibC(1,i1,k) = (1+ abs((ZAR(i1,1)-CD.W)/(ZAR(i1,1)+CD.W)))/(1- abs((ZAR(i1,1)-CD.W)/(ZAR(i1,1)+CD.W)));
    Ksv.VibD(1,i1,k) = (1+ abs((ZAR(i1,2)-CD.W)/(ZAR(i1,2)+CD.W)))/(1- abs((ZAR(i1,2)-CD.W)/(ZAR(i1,2)+CD.W)));

end
    %Сопротивления в точке C и D
    ZBx.C(k) = 1/((1/(Z.C(1)))  + (1/(Z.C(2)))  + (1/(Z.C(3))));
    ZBx.D(k)=  1/((1/(Z.D(1)))  + (1/(Z.D(2)))  + (1/(Z.D(3))));
    %%КСВ
    Ksv.C(k) = (1+ abs((ZBx.C(k)-B.W)/(ZBx.C(k)+B.W)))/(1- abs((ZBx.C(k)-B.W)/(ZBx.C(k)+B.W)));
    Ksv.D(k) = (1+ abs((ZBx.D(k)-B.W)/(ZBx.D(k)+B.W)))/(1- abs((ZBx.D(k)-B.W)/(ZBx.D(k)+B.W)));
    
       %%
                                 %Блок B 
    %Сопртоивления с учетом линии
    Z.B(1) = (B.W* (ZBx.C(k)+j*B.W*tan(B.beta_coax(k)*B.L))/(B.W+j*ZBx.C(k)*tan(B.beta_coax(k)*B.L)));
    Z.B(3) = Z.B(1);
    Z.B(2) = (B.W* (ZBx.D(k)+j*B.W*tan(B.beta_coax(k)*B.L))/(B.W+j*ZBx.D(k)*tan(B.beta_coax(k)*B.L)));
    %Сопротивления в точке B
    ZBx.B(k) = 1/(  (1/(Z.B(1)))  + (1/(Z.B(2)))  + (1/(Z.B(3))));

    % влияние четвертьволнового трансформатора  
    ZBx.B(k)=B.Wtr *(ZBx.B(k) + j*B.Wtr*tan(6.28*B.Ltr/B.lamB(k)))/(B.Wtr+j*ZBx.B(k)*tan(6.28*B.Ltr/B.lamB(k)));  
%     %Влияние трехступенчатого трансформатора сопротивления                               
%     B.Z3= B.W3*(ZBx.B(k)+j*B.W3*tan(Freq.beta_coax(k)*Freq.L))/(B.W3+j*ZBx.B(k)*tan(Freq.beta_coax(k)*Freq.L));
%     B.Z2= B.W2*(B.Z3+j*B.W2*tan(Freq.beta_coax(k)*Freq.L))/(B.W2+j*B.Z3*tan(Freq.beta_coax(k)*Freq.L));
%     ZBx.B(k) =B.W1*(B.Z2+j*B.W1*tan(Freq.beta_coax(k)*Freq.L))/(B.W1+j*B.Z2*tan(Freq.beta_coax(k)*Freq.L));   
    %%КСВ 
    Ksv.B(k) = (1+ abs((ZBx.B(k)-A.W)/(ZBx.B(k)+A.W)))/(1- abs((ZBx.B(k)-A.W)/(ZBx.B(k)+A.W)));
 
 %%
                                %Блок A
    %Сопртоивления с учетом линии
    Z.A(1)= (A.W* (ZBx.B(k)+j*A.W*tan(A.beta_coax(k)*A.L))/(A.W+j*ZBx.B(k)*tan(A.beta_coax(k)*A.L)));
    %Сопротивления в точке A
    ZBx.A(k)= 1/ (2/Z.A(1));  
    
%     Влияние 1/4 трансформатора 
ZBx.A(k) = A.Wtr* (ZBx.A(k) + j*A.Wtr*tan(6.28*A.Ltr/B.lamB(k)))/(A.Wtr + j*ZBx.A(k)*tan(6.28*A.Ltr/B.lamB(k)));

    
%     %%Влияние трансформатора сопротивления  
%     A.Z3= A.W3*(ZBx.A(k)+j*A.W3*tan(Freq.beta_coax(k)*Freq.L))/(A.W3+j*ZBx.A(k)*tan(Freq.beta_coax(k)*Freq.L));
%     A.Z2= A.W2*(A.Z3+j*A.W2*tan(Freq.beta_coax(k)*Freq.L))/(A.W2+j*A.Z3*tan(Freq.beta_coax(k)*Freq.L));                               
%     %Сопротивления в точке A
%     ZBx.A(k) =A.W1*(A.Z2+j*A.W1*tan(Freq.beta_coax(k)*Freq.L))/(A.W1+j*A.Z2*tan(Freq.beta_coax(k)*Freq.L));   
    %%КСВ
    Ksv.A(k) = (1+ abs((ZBx.A(k)-A.W)/(ZBx.A(k)+A.W)))/(1- abs((ZBx.A(k)-A.W)/(ZBx.A(k)+A.W)));

 %%
   
    
    
ZORRO(k) = sum(sum(ZAR));
     if k==length(Freq.lambda0)
ZORROZAR = ZAR;
     end
     clear Z ZA ZAR ZORR Freq.beta_coax(k) Freq.L;
   
disp(k);
    



 end


for i = 1:length(Ksv.VibC) 
Ksv.V1(i)=Ksv.VibC(1,1,i);
Ksv.V4(i)=Ksv.VibC(1,2,i);
Ksv.V7(i)=Ksv.VibC(1,3,i);
Ksv.V2(i)=Ksv.VibD(1,1,i);
Ksv.V5(i)=Ksv.VibD(1,2,i);
Ksv.V8(i)=Ksv.VibD(1,3,i);
end    


 

%% Figures
createfigure(Freq.F,Ksv.C);
createfigure(Freq.F,Ksv.D);
createfigure(Freq.F,Ksv.B);
createfigure(Freq.F,Ksv.A);
figure,plot(Freq.F,Ksv.A,Freq.F,Ksv.B,Freq.F,Ksv.D,Freq.F,Ksv.C)

figure,plot(Freq.F,Ksv.V1,Freq.F,Ksv.V4,Freq.F,Ksv.V7);
figure,plot(Freq.F,Ksv.V2,Freq.F,Ksv.V5,Freq.F,Ksv.V8);
