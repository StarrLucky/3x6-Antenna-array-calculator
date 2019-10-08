function createaxes(Parent1, X1, Y1)
%CREATEAXES1(PARENT1,X1,Y1)
%  PARENT1:  axes parent
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 01-Feb-2013 10:37:31

% Create axes
axes1 = axes('Parent',Parent1,'YGrid','on',...
    'XTick',[0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400],...
    'XMinorTick','on',...
    'XGrid','on');
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1,'LineWidth',1.5);

% Create xlabel
xlabel({'\theta'});

% Create ylabel
ylabel({'F(\theta)'});
