% close all;
%%

Freq.lambda0 = 0.19:-0.001:0.175;        %Предполагаемый диапазон длин работы антенны

Freq.lambda0_= 0.18;        %Рабочая длина волны
kk = 1.19;                  %коэффициент укорочения


%%
% Freq.lambda0=0.18;
W=50;  
Freq.C=299792458;  
Freq.F = Freq.C./Freq.lambda0; 
Freq.lam1 = Freq.lambda0./kk; %Длина волны в коакс. фидере
Freq.lam1_0 =Freq.lambda0_/kk;
%%
CD.L = 0.08; B.L = 0; A.L = 0;          %Длины кабеля в секциях
%Волновые споротивления трехступенчатых трансформаторов
B.W1 = 43.63; B.W2 = 33.23; B.W3 = 25.31;    
A.W1 = 44.54; A.W2 = 35.35; A.W3 = 28.06;
Freq.L = 0.1513/4;    % Длина секции трехступенчатого трансформатора




%%

for k=1:length(Freq.lambda0) 

   
    % k = 1;

   
                                              Freq.beta_coax(k) =6.28./Freq.lam1(k);
   %%
                                              lam = Freq.lambda0(k);   %#ok<NASGU> %Global

                                              

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
    %Сопротивления с учетом линии
    Z.C(i1) = W* (ZAR(i1,1)+j*W*tan(Freq.beta_coax(k)*CD.L))/(W+j*ZAR(i1,1)*tan(Freq.beta_coax(k)*CD.L));  %Трансформация сопротивления в линии
    Z.C(i1) =1/ ((1/(Z.C(i1))) + (1/ (j*W*tan(Freq.beta_coax(k)*(Freq.lam1_0/4)))));  %Влияние четвертьволнового стакана
   
    Z.D(i1) = (W* (ZAR(i1,2)+j*W*tan(Freq.beta_coax(k)*CD.L))/(W+j*ZAR(i1,2)*tan(Freq.beta_coax(k)*CD.L)));
    Z.D(i1) = 1/ ((1/(Z.D(i1))) + (1/ (j*W*tan(Freq.beta_coax(k)*(Freq.lam1_0/4)))));

end
    %Сопротивления в точке C и D
    ZBx.C(k) = 1/((1/(Z.C(1)))  + (1/(Z.C(2)))  + (1/(Z.C(3))));
    ZBx.D(k)=  1/((1/(Z.D(1)))  + (1/(Z.D(2)))  + (1/(Z.D(3))));
    %%КСВ
    Ksv.C(k) = (1+ abs((ZBx.C(k)-W)/(ZBx.C(k)+W)))/(1- abs((ZBx.C(k)-W)/(ZBx.C(k)+W)));
    Ksv.D(k) = (1+ abs((ZBx.D(k)-W)/(ZBx.D(k)+W)))/(1- abs((ZBx.D(k)-W)/(ZBx.D(k)+W)));
 
       %%
                                 %Блок B 
    %Сопртоивления с учетом линии
    Z.B(1) = (W* (ZBx.C(k)+j*W*tan(Freq.beta_coax(k)*B.L))/(W+j*ZBx.C(k)*tan(Freq.beta_coax(k)*B.L)));
    Z.B(3) = Z.B(1);
    Z.B(2) = (W* (ZBx.D(k)+j*W*tan(Freq.beta_coax(k)*B.L))/(W+j*ZBx.D(k)*tan(Freq.beta_coax(k)*B.L)));
    %Сопротивления в точке B
    ZBx.B(k) = 1/(  (1/(Z.B(1)))  + (1/(Z.B(2)))  + (1/(Z.B(3))));
    %%Влияние трансформатора сопротивления                               
    B.Z3= B.W3*(ZBx.B(k)+j*B.W3*tan(Freq.beta_coax(k)*Freq.L))/(B.W3+j*ZBx.B(k)*tan(Freq.beta_coax(k)*Freq.L));
    B.Z2= B.W2*(B.Z3+j*B.W2*tan(Freq.beta_coax(k)*Freq.L))/(B.W2+j*B.Z3*tan(Freq.beta_coax(k)*Freq.L));
    %Сопротивления в точке B
    ZBx.B(k) =B.W1*(B.Z2+j*B.W1*tan(Freq.beta_coax(k)*Freq.L))/(B.W1+j*B.Z2*tan(Freq.beta_coax(k)*Freq.L));   
    %%КСВ
    Ksv.B(k) = (1+ abs((ZBx.B(k)-W)/(ZBx.B(k)+W)))/(1- abs((ZBx.B(k)-W)/(ZBx.B(k)+W)));
 
 %%
                                %Блок A
    %Сопртоивления с учетом линии
    Z.A(1)= (W* (ZBx.B(k)+j*W*tan(Freq.beta_coax(k)*A.L))/(W+j*ZBx.B(k)*tan(Freq.beta_coax(k)*A.L)));
    %Сопротивления в точке A
    ZBx.A(k)= 1/ (2/Z.A(1));                                                                                          
    %%Влияние трансформатора сопротивления  
    A.Z3= A.W3*(ZBx.A(k)+j*A.W3*tan(Freq.beta_coax(k)*Freq.L))/(A.W3+j*ZBx.A(k)*tan(Freq.beta_coax(k)*Freq.L));
    A.Z2= A.W2*(A.Z3+j*A.W2*tan(Freq.beta_coax(k)*Freq.L))/(A.W2+j*A.Z3*tan(Freq.beta_coax(k)*Freq.L));                               
    %Сопротивления в точке A
    ZBx.A(k) =A.W1*(A.Z2+j*A.W1*tan(Freq.beta_coax(k)*Freq.L))/(A.W1+j*A.Z2*tan(Freq.beta_coax(k)*Freq.L));   
    %%КСВ
    Ksv.A(k) = (1+ abs((ZBx.A(k)-W)/(ZBx.A(k)+W)))/(1- abs((ZBx.A(k)-W)/(ZBx.A(k)+W)));

 %%
   
    
    
ZORRO(k) = sum(sum(ZAR));
     if k==length(Freq.lambda0)
ZORROZAR = ZAR;
     end
     clear Z ZA ZAR ZORR Freq.beta_coax(k) Freq.L;
   
disp(k);
    
end



%% Figures
createfigure(Freq.F,Ksv.C);
createfigure(Freq.F,Ksv.D);
createfigure(Freq.F,Ksv.B);
createfigure(Freq.F,Ksv.A);
figure,plot(Freq.F,Ksv.A,Freq.F,Ksv.B,Freq.F,Ksv.D,Freq.F,Ksv.C)