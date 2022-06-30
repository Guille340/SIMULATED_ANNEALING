% DESCRIPTION
% The camelback function has several local minima plus two global minima
% at [-0.0898, 0.7126] and [0.0898, -0.7126], with an associated value (cost)
% of -1.0316. The following example attempts to find these solutions using
% the simulated annealing function SIMANNEAL.m

% Options Structure
Options = simannealOptions('InitialTemperature',100,'CoolingFactor',0.9,...
    'minTemperature',0.01,'MaxAccept',20,'plot','on');

% Anneal
camelFun = @(x,y) (4 - 2.1*x.^2 + x.^4/3).*x.^2 + x.*y + 4*(y.^2 - 1).*y.^2;
lossFun = @(p) camelFun(p(1),p(2));
[xi,Ci,History] = simanneal(lossFun,[5 -5],[-10 -10],[10 10],Options);

% Plot Results
figure
hold on
[x,y] = meshgrid(-10:0.1:10,-10:0.1:10);
z = (4 - 2.1*x.^2 + x.^4/3).*x.^2 + x.*y + 4*(-1 + y.^2).*y.^2;
surf(x,y,z)
hax = gca;
hax.Children.EdgeColor = 'none';
plot3(History.solutions(1,:),History.solutions(2,:),History.costs,'ko',...
    'MarkerSize',2,'MarkerFaceColor','k','Color','none')

