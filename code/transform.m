% % Расчет четвертьволнового трансформатора
%%
% tr_Rn =24.6334 + 4.2659i;
% tr_lam = 0.1513;  %!
% tr_W=50;
% tr_L = tr_lam/4;
% tr_b = 2*pi/tr_lam;
% 
% %%
% tr_W2 = sqrt(tr_Rn*tr_W); tr_W3=(tr_Rn^2 *tr_W2)^(1/3); 
% tr_W1 = tr_Rn*tr_W/tr_W3;
% 
% tr_W2=abs(tr_W2);
% tr_W1=abs(tr_W1);tr_W3=abs(tr_W3);
% tr_Z3= tr_W3*(tr_Rn+j*tr_W3*tan(tr_b*tr_L))/(tr_W3+j*tr_Rn*tan(tr_b*tr_L));
% tr_Z2= tr_W2*(tr_Z3+j*tr_W2*tan(tr_b*tr_L))/(tr_W2+j*tr_Z3*tan(tr_b*tr_L));
% tr_Z1=tr_W1*(tr_Z2+j*tr_W1*tan(tr_b*tr_L))/(tr_W1+j*tr_Z2*tan(tr_b*tr_L));
% tr_KSV = (1+ abs((tr_Z1-tr_W)/(tr_Z1+tr_W)))/(1- abs((tr_Z1-tr_W)/(tr_Z1+tr_W)));


% W = 50;
% Z = 28.3913118624425 - 9.90650543115715i;
% (1+ abs((Z-W)/(Z+W)))/(1- abs((Z-W)/(Z+W))); %#ok<VUNUS>
%  


