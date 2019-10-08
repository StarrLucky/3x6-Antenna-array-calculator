                                                close all;
%                                                 clear all;
%% Диаграммы направленности
%%
                                        x3d = 0;    %Диаграммы в 3D
                                        x2d =1;     % Диаграммы в 2D
lam = 0.18;       %Рабочая
% lam =0.1968;    %Нижняя
% lam =0.1637;    %Верхняя
%%
                                          teta = -pi/2:0.01:pi/2; 
%                                           teta = 0.6109:0.0001:0.7505; 
%0.3491
%%
k = 2*pi/lam;
Nx = 3;  %
Ny = 6;  % 
dx =0.09;   %(Dx) расстояние между соседними излучателями по оси OX
dy = 0.045; % (Dy) расстояние между соседними излучателями по оси OX
dp = 0.045;
l=0.045;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ДН антенной решетки в E и H плоскостях


%E.Fp = (pi/4 * (cos(teta) - 1));
E.Fp = 2*sin(k*dp*cos(teta));
phi=0;
E.F0 = ( cos((pi/2) * sin(teta))./ cos(teta));
% E.F0 = (cos(k*l*cos(teta+pi/2))- cos(k*l))./(1-cos(k*l)*sin(teta+pi/2));
E.Fa = abs(sin(Nx*(k*dx*sin(teta)*cos(phi))/2)./(sin(0.5*(k*dx*sin(teta)*cos(phi)))));
%c
phi=pi/2;
H.F0=linspace(1,1,length(teta));
H.Fa = abs(sin(Ny*(k*dy*sin(teta)*sin(phi))/2)./(sin(0.5*(k*dy*sin(teta)*sin(phi)))));
H.Fz =1;

%%
E.F = (E.Fp.*abs(E.F0).*abs(E.Fa));
H.F = (E.Fp.*abs(H.Fa)); 
      
% for i = 1:length(E.F)
%    if E.F(i) < 0 
%    E.F(i) = NaN;
%    end  
% end
% for i = 1:length(H.F)
%    if H.F(i) < 0 
%    H.F(i) = NaN;
%    end  
% end
   
%%   
E.max = max(abs(E.F));
H.max = max(abs(H.F));



%%
if x3d ==1 
%%3d plotting 
Surf.Theta = -pi/2:0.05:pi/2;
Surf.Phi= -pi:0.1:pi; 

Surf.max = 0;Surf.min = 0; Surf.tetamax = 0;
Surf.phimax = 0; Surf.amp=0;
                     %% %% Rho, Phi, Theta --> x,y,z 
                     Surf.S = 0;
                     Surf.percent = length(Surf.Phi)*length(Surf.Theta)/100;
for k0=1:length(Surf.Phi)
for k1=1:length(Surf.Theta) 
Surf.S=Surf.S+1;
%Множитель E+H
Surf.Tempo =abs(sin(Nx*(k*dx*sin(Surf.Theta(k1))*cos(Surf.Phi(k0)))/2)/(sin(0.5*(k*dx*sin(Surf.Theta(k1))*cos(Surf.Phi(k0))))))*(sin(Ny*(k*dy*sin(Surf.Theta(k1))*sin(Surf.Phi(k0)))/2)/(sin(0.5*(k*dy*sin(Surf.Theta(k1))*sin(Surf.Phi(k0))))));



% Surf.Tempo =(sin(Nx*(k*dx*sin(Surf.Theta(k1))*cos(Surf.Phi(k0)))/2)/(sin(0.5*(k*dx*sin(Surf.Theta(k1))*cos(Surf.Phi(k0))))));
% Surf.Tempo =(sin(Ny*(k*dy*sin(Surf.Theta(k1))*sin(Surf.Phi(k0)))/2)/(sin(0.5*(k*dy*sin(Surf.Theta(k1))*sin(Surf.Phi(k0))))));

%Вибратор в E
Surf.Tempo= Surf.Tempo* (abs( cos((pi/2) * sin(Surf.Theta(k1)))./ cos(Surf.Theta(k1))));  

%Вибратор в H
% Surf.Tempo= Surf.Tempo* 1;

Surf.Tempo= Surf.Tempo* (2 * (sin(k*dp*cos(Surf.Theta(k1)))));  
%Земля
% Surf.Tempo = Surf.Tempo* 4*sin(k*sin(Surf.Theta(k1))*Hsr);

Surf.amp=Surf.Tempo;

Surf.Tempo= abs(Surf.Tempo);
   
    Surf.x(Surf.S) = (Surf.Tempo) * sin(Surf.Theta(k1))*cos(Surf.Phi(k0));
    Surf.y(Surf.S) = (Surf.Tempo) * sin(Surf.Theta(k1))*sin(Surf.Phi(k0));
    Surf.z(Surf.S) = (Surf.Tempo) * cos(Surf.Theta(k1));

    
    if Surf.amp > Surf.max
        Surf.max=Surf.amp;
        Surf.tetamax = Surf.Theta;
        Surf.phimax =Surf.Phi;
    end
     if Surf.amp < Surf.min
        Surf.min=Surf.amp;
     end
    
    clear Surf.Tempo Surf.amp;
    
end
disp(Surf.S/Surf.percent); 
end

end




%% Графики
if x2d == 1 
figure, polar(teta, E.Fa);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off
figure, polar(teta, E.F/E.max);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off

hold on
createaxes(figure,radtodeg((teta)),(E.F./E.max));
hold off

hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off

figure, polar(teta, H.Fa);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off

figure, polar(teta, H.F/H.max);
hold on
xlabel({'\theta'});
ylabel({'F(\theta)'});
hold off



hold all
createaxes(figure,radtodeg((teta)),(H.F./H.max));
hold off
%

% plot3();
tetadeg = radtodeg(teta);
% teta3=teta-(pi);
H.Fnorm = abs(H.F/H.max);
E.Fnorm = abs(E.F/E.max);
end


if x3d==1 
    figure,plot3( Surf.x, Surf.y, Surf.z);
end
