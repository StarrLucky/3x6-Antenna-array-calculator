% % close all
% % clear all
% % 
% % 
% % % y = linspace(1,10);
% % % 
% % % s= sin(y);
% % % s(10) = 5;
% % % 
% % % 
% % % figure,plot(s);
% % % 
% % % %figure,plot(s./y);
% % % 
% % % 
% % % 
% % % l = 0;
% % % 
% % % r = length(y);
% % % h(1:r) = 0;
% % % g(1:1)=0;
% % % for i1 =1:r
% % %     
% % %   l = (r+1) - i1;
% % %   
% % %   h(i1) = s(i1)/y(l)^2;
% % %   g(i1) = s(i1)/y(i1)^2;
% % %     
% % % end
% % % 
% % % figure,plot(h);
% % % figure,plot(g);
% % % 
% % % 
% % % a = g-s;
% % % figure,plot(a);
% % 
% % 
% % % hex2dec('0x07');
% % % dec2hex(1024);
% % 
% % 
% % clear all;
% % 
% % % 
% % % 
% % % E=2.25;
% % % % W = 138* (sqrt(1/E))*(log10(0.0048/0.00137));
% % % % 
% % % 
% % % syms d;
% % % dsolve ('log10(d) - 2.84 =0')
% % % 
% % % % 
% % % D = 12.9;
% % % 
% % % (D+1)/(64.8^2  + 1);
% %  syms d;
% % fsolve('276*log10( (12.9/d) + sqrt( (12.9/d) -1))=500',d)
% % 
% % k=2*pi/0.18;
% % l=0.0445;
% % Rs =81.0689;
% % Rsp = 81.0442;
% % W=333.4304;
% % 
% % ZZ= (Rsp/(((Rsp/W)^2)+(sin(k*l))^2)) - (j*(W/2) * ((sin(2*k*l))^2)/((Rsp/W)^2 - (sin(k*l))^2));
% 
% % Wa = 333.4304;
% % 
% % k = 2* pi / 0.18;
% % l = 0.0413;
% % 
% % 
% % Rsp = 81.0689 * (sin(2*3.14*0.0413/0.18))^2;
% % 
% % Za = (Rsp/( (Rsp/Wa)^2 + (sin(k*l))^2)) - j*(Wa/2)*(sin(2*k*l)/( (Rsp/Wa)^2 + (sin(k*l))^2));
% 
% 
% 148.972933955072 + 33.9374310356292i
% % 148.972933955072 + 33.9374310356292i
% % Part1
% 
ZN = (162.016966621278 + 76.7915467324876i);
W= 150;
lam = 0.1513;


l0=0.0205 + 0.5;


l0 = 0.0948+0.5;

l = l0 *lam;


 
Zn = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));

% 
% % 
for i1 = 1:3
ZN = ZAR(i1,1);
OPT(i1,1) = W*(ZN+W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
Ksvvib(i1,1) =  (1+ abs((ZN-W)/(ZN+W)))/(1- abs((ZN-W)/(ZN+W)));
 ZN = ZAR(i1,2);
OPT(i1,2) = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
Ksvvib(i1,2) =  (1+ abs((ZN-W)/(ZN+W)))/(1- abs((ZN-W)/(ZN+W)));

end
% % Rsum
OPT(1,3) = 1/( (1/OPT(1,1))+(1/(OPT(2,1)))+ (1/(OPT(3,1))));
OPT( 2,3) = 1/((1/OPT(1,2))+(1/(OPT(2,2)))+ (1/(OPT(3,2))));
% % SWR

W=75;
  OPT(1,4) =  (1+ abs((OPT(1,3)-W)/(OPT(1,3)+W)))/(1- abs((OPT(1,3)-W)/(OPT(1,3)+W)));
  
OPT(2,4)=(1+ abs((OPT(2,3)-W)/(OPT(2,3)+W)))/(1- abs((OPT(2,3)-W)/(OPT(2,3)+W))); 




% Part2
OPT(3,3) = OPT(1,3);
W= 75;
lam = 0.1184;

l0=0.0830;

l0=0.5068;


l = l0*lam;

for i1 = 1:3
ZN = OPT(i1,3);
OPT(i1,5) = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
end
% % % % Rsum
OPT(1,6) = 1/( (1/OPT(1,5))+(1/(OPT(2,5)))+ (1/(OPT(3,5))));

% % % % SWR
W=50;

  OPT(1,7) =  (1+ abs((OPT(1,6)-W)/(OPT(1,6)+W)))/(1- abs((OPT(1,6)-W)/(OPT(1,6)+W)));
%   

% %%Part3
% 
% OPT(2,6) = OPT(1,6);
% 
% ZN = Z1;
% 
% W= 50;
% lam = 0.1513;
%  llam = 0.1340 + 0.247;
% l = lam *llam;
% % l= 0.06456;
% % l= 0.064;
% 
% 
% Zn = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
% for i1 = 1:3
% ZN = OPT(i1,3);
% OPT(i1,5) = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
% end
% %Rsum
% OPT(1,6) = 1/( (1/OPT(1,5))+(1/(OPT(2,5)))+ (1/(OPT(3,5))));
% 
% % SWR
%   OPT(1,7) =  (1+ abs((OPT(1,6)-W)/(OPT(1,6)+W)))/(1- abs((OPT(1,6)-W)/(OPT(1,6)+W)));
% %   
% % 
% % 
% % 
% % %
% % % %part 5 блок A
% % 
% % % transform;
% % 
% % OPT(1,8) = Z1;OPT(2,8) = Z1;
% % 
% % 
% % % ZN = Z1;
% % 
% % W= 50;
% % lam = 0.12;
% %  llam = 0.292;
% % l = lam *llam;
% % % l= 0.06456;
% % % l= 0.064;
% % 
% % 
% % Zn = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));
% % for i1 = 1:2
% % ZN = OPT(i1,8);
% % OPT(i1,9) = W*(ZN +W*j*tan(6.28*l/lam))/(W+ZN*j*tan(6.28*l/lam));  %Сопротивление с линией
% % end
% % %Rsum
% % OPT(1,10) = 1/( (1/OPT(1,9))+(1/(OPT(2,9))));
% % 
% % % SWR
% %   OPT(1,11) =  (1+ abs((OPT(1,10)-W)/(OPT(1,10)+W)))/(1- abs((OPT(1,10)-W)/(OPT(1,10)+W)));
% %   
% 
% 
% % (ZN +W*j*tan(6.28*l/lam))
% % (W+ZN*j*tan(6.28*l/lam))


% result = dblquad(@tests,0,pi,0,2*pi)


% 
% function  out = tests(teta,phi)
% 
% out = (sin(Ny*(k*dx*sin(teta)*sin(phi))/2)./(sin(0.5*(k*dy*sin(teta)*sin(phi)))))*(sin(Nx*(k*dx*sin(teta)*cos(phi))/2)./(sin(0.5*(k*dx*sin(teta)*cos(phi))))); 
% 
% end
% 






% 
% function out = tests(teta,phi)
% out = y*sin(x) + x*cos(y); 
% Ny = 6; Nx=3; dx=0.09; dy=0.045;
% k = 37.17;
%  out = (sin(Ny*(k*dx*sin(teta)*sin(phi))/2)/(sin(0.5*(k*dy*sin(teta)*sin(phi)))))*(sin(Nx*(k*dx*sin(teta)*cos(phi))/2)/(sin(0.5*(k*dx*sin(teta)*cos(phi))))) *sin(teta); 
% 
% end
% 
% 
% syms a
% % solve('(1/(1+( (0.18/2*a)*(1/(0.0413*(2*3.14*0.0063/a))))^2))=0.07',a)
%  
% 
% 
% 
% 
% W=333.43;
% l=0.045;
% k=2*pi/0.18;
% ZIZAR(6,3) = 0; 
% i = sqrt(-1);
% for i1 = 1:6
%     for j1 = 1:3
%  Z = ZAR(i1,j1);       
% ZIZAR(i1,j1)=-(1/2)*(-W+sqrt(W^2-4*Z*i*W*sin(k*l)*cos(k*l)-4*Z^2*sin(k*l)^2))*W/Z;
% clear Z;
%     end
% end
% 