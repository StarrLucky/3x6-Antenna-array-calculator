function createfigure(X1, Y1)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 04-Feb-2013 09:12:13

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YTick',[1 1.5 2 2.5 3 3.5 4],...
    'YMinorTick','on',...
    'YGrid','on',...
    'XTick',[1500000000 1540000000 1580000000 1620000000 1660000000 1700000000 1740000000 1780000000 1820000000 1850000000],...
    'XMinorTick','on',...
    'XGrid','on');
box(axes1,'on');
hold(axes1,'all');

% Create plot
plot(X1,Y1);

% Create xlabel
xlabel({'F, ��'});

% Create ylabel
ylabel({'���'});
