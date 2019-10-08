%%Dn_v_over10
close all;
clear all;

x2d =1;     % Диаграммы в 2D
%ДН антенной решетки с полуволновыми вибраторами
%%
% lam = 0.18;       %Рабочая
% lam =0.189;    %Нижняя
lam =0.1756;    %Верхняя
% 
k = 2*pi/lam;
Nx = 3;  %
Ny = 6;  % 
dx =0.09;   %(Dx) расстояние между соседними излучателями по оси OX
dy = 0.045; % (Dy) расстояние между соседними излучателями по оси OX
dr = 0.045;
l = lam/4;


theta = 0:0.01:2*pi;
delta = 0:0.01:2*pi;

E.F0 = (cos(k*l*cos(theta))-cos(k*l))./(1-cos(k*l)*sin(theta));
E.Fp = 2*sin(k*dr*cos(theta-(pi/2)));
E.Fa = sin((1/2)*Nx*k*dx*cos(theta))./(sin((1/2)*k*dx*cos(theta)));

E.F = E.Fa.*E.F0.*E.Fp;
E.max = max(abs(E.F));

% 
% for i = 1:length(E.F)
%    if E.F(i) < 0 
%    E.F(i) = NaN;
%    end  
% end


H.Fp = 2*sin(k*dr*cos(delta-(1/2)*pi));
H.Fa = sin((1/2)*Ny*k*dy*cos(delta))./(sin((1/2)*k*dy*cos(delta)));
H.F =( H.Fp.*H.Fa);

H.max = max(abs(H.F));


% for i = 1:length(H.F)
%    if H.F(i) < 0 
%    H.F(i) = NaN;
%    end  
% end



thetadeg=radtodeg(theta);
deltadeg= radtodeg(delta);



% 
% figure,polar(theta,Fp);
% figure,polar(theta,Efc);
% figure,polar(theta,FE./maxE);
% 
% 
% figure,polar(delta,HFref);
% figure,polar(delta,Hfc);
% figure,polar(delta,FH./maxH);





%% Графики
if x2d == 1 
figure, polar(theta, E.Fa);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off
figure, polar(theta, E.F/E.max);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off

hold on
createaxes(figure,radtodeg(theta),(E.F./E.max));
hold off

hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off

figure, polar(delta, H.Fa);
hold on
xlabel({'\Delta'});
ylabel({'F(\Delta)'});
hold off

figure, polar(delta, H.F/H.max);
hold on
xlabel({'\Delta'});
ylabel({'F(\Delta)'});
hold off



hold all
createaxes(figure,radtodeg(delta),(H.F./H.max));
hold off
%

% plot3();
thetadeg = radtodeg(theta);
% teta3=teta-(pi);
H.Fnorm = (H.F/H.max);
E.Fnorm = (E.F/E.max);
end


