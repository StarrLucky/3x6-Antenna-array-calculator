%%              ¬заимные сопротивлени€ вибраторов решетки

%%
ZAR(m,n) = 0;           %Ќаведенные сопротивлени€ вибраторов решетки
ZA(m,n,m*n*2) = 0;      % ћассив сопротивлений R,jX. ...
                        %(i,j, 1..m*n) - взаимные сопротивлени€ i,j вибраторов антенной решетки
                        %(i,j, m*n..2*m*n) - взаимные сопротивлени€ i,j вибратора и вибраторов
                        % рефлектора
%%


for ii =1:m         % 1:(round(m/2))  
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
%%—уммирование наведенных сопротивлений вибратора от других вибраторов
%%и рефлектора
 ZAR(ii,jj) = ZAR(ii,jj)+ ZA(ii,jj,S) +  (ZA(ii,jj,S+m*n)*(-1));

end
end
end
end

clear ii i0 jj j0;
clear hh bb dd;