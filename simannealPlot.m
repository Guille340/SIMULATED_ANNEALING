%  simannealPlot(History)
%
%  DESCRIPTION
%  Plots the results from SIMANNEAL.m returned as a HISTORY structure. The 
%  function is run is the 'Plot' option in the OPTIONS structure is enabled
%  (OPTIONS.PLOT = 'on').
%
%  There are four types of plot:
%  - "Solution History": evolution of the value of each individual variable
%    with time.
%  - "Cost + Temperature History": evolution of the cost and temperature
%    with time. The values are scaled with respect to the left y-axis and the
%    right y-axis.
%  - "Cost Difference History": evolution of the cost difference with time.
%  - "Cost per Solution": cost as a function of the value of each individual
%    variable.
%
%  INPUT ARGUMENTS
%  - History: structure containing the optimisation parameters for every 
%    iteration happening in the annealing process. Generated with function
%    SIMANNEAL.m.
%
%  OUTPUT ARGUMENTS
%  - None
%
%  FUNCTION CALL
%  1. simannealPlot(History)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  See also SIMANNEALOPTIONS, SIMANNEAL

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Jun 2022

function simannealPlot(History,varargin)

narginchk(1,2)

% Input Arguments
if nargin == 1
    figDir = fullfile(pwd,'Figures');
else % nargin == 2;
    figDir = varargin{1};
end

% Error Control ('History')
if ~isstruct(History)
    error('Input argument HISTORY must be a valid structure')
end

% Error Control ('figDir')
if ~ischar(figDir) || ~ismember('\',figDir)
    error('VARARGIN{1} is not a valid character vector for a folder')
end

% Identify Success, Fail and Rejection Solutions
iAccept = History.accepts == 1;
iSuccess = History.accepts == -1;
iReject = History.accepts == 0;
nIter = length(History.costs); % iteration vector
iter = 1:nIter;

% Create 'Figures' Directory

mkdir(figDir)

% Solution History
nVariables = size(History.solutions,1);
for n = 1:nVariables
    figure
    hold on
    h(1) = scatter(iter(iReject),History.solutions(n,iReject),20,...
        [0.85 0.85 0.85],'filled');
    h(2) = scatter(iter(iSuccess),History.solutions(n,iSuccess),20,...
        [0.5 0.5 0.5],'filled');
    h(3) = scatter(iter(iAccept),History.solutions(n,iAccept),20,[0 0 0],'filled');
    title(sprintf(['\\bfSolution History of Variable %d \\rm'...
        '(Best Solution = %0.3g)'],n,History.bestSolution(n)))
    xlabel('Iteration','FontSize',10)
    ylabel(sprintf('Value of Variable x(%d)',n),'FontSize',10)
    pbaspect([1 1 1])
    legend(h,{'rejected','accepted ("bad")','accepted ("good")',},...
        'Location','SouthEast','FontSize',10)
%     ymin = History.bestSolution(n) - 2*std(History.solutions(n,:));
%     ymax = History.bestSolution(n) + 2*std(History.solutions(n,:));
    ymin = min(History.solutions(n,:));
    ymax = max(History.solutions(n,:));
    if ymin ~= ymax, ylim([ymin ymax]); end
    xlim([1 nIter])
    clear h 
    box on
    set(gcf,'Color',[1 1 1])
    set(gcf,'PaperPositionMode','auto','InvertHardcopy','off')
    set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.8])
    figPath = fullfile(figDir,sprintf('Solution History (Var%d)',n));
    print(strcat(figPath,'.png'),'-dpng','-r250');
    savefig(strcat(figPath,'.fig'))
    close(gcf)
end

% Cost + Temperature History
figure
hold on

yyaxis right
h = plot(iter,History.temperatures,'k','LineWidth',1.5);
ylim([min(History.temperatures) max(History.temperatures)])
ylabel('Temperature','FontSize',10)

yyaxis left
h(1) = scatter(iter(iReject),History.costs(iReject),20,[0.85 0.85 0.85],'filled');
h(2) = scatter(iter(iSuccess),History.costs(iSuccess),20,[0.5 0.5 0.5],'filled');
h(3) = scatter(iter(iAccept),History.costs(iAccept),20,[0 0 0],'filled');

% ymin = History.bestCost - 0.1*std(History.costs);
% ymax = History.bestCost + 0.5*std(History.costs);
ymin = min(History.costs);
ymax = max(History.costs);
if ymin ~= ymax, ylim([ymin ymax]); end
xlim([1 nIter])
ylabel('Cost C','FontSize',10)
pbaspect([1 1 1])
legend(h,{'rejected','accepted ("bad")','accepted ("good")'},...
    'Location','NorthEast','FontSize',10)
clear h    
title('Cost and Temperature History')
xlabel('Iteration','FontSize',10)
hax = gca;
hax.YAxis(1).Color = [0 0 0];
hax.YAxis(2).Color = [0 0 0];
box on  
set(gcf,'Color',[1 1 1])
set(gcf,'PaperPositionMode','auto','InvertHardcopy','off')
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.8])
figPath = fullfile(figDir,'Cost & Temperature History');
print(strcat(figPath,'.png'),'-dpng','-r250');
savefig(strcat(figPath,'.fig'))
close(gcf)

% Cost Difference History
figure
hold on
h(1) = scatter(iter(iReject),History.costDifferences(iReject),20,[0.85 0.85 0.85],'filled');
h(2) = scatter(iter(iSuccess),History.costDifferences(iSuccess),20,[0.5 0.5 0.5],'filled');
h(3) = scatter(iter(iAccept),History.costDifferences(iAccept),20,[0 0 0],'filled');
ylim([min(History.costDifferences) max(History.costDifferences)])
title('Cost Difference History')
xlabel('Iteration','FontSize',10)
ylabel('Cost Difference \DeltaC','FontSize',10)
pbaspect([1 1 1])
legend(h,{'rejected','accepted ("bad")','accepted ("good")'},...
    'Location','SouthEast','FontSize',10)
% ymin = -0.2*std(History.costDifferences);
% ymax = 0.2*std(History.costDifferences);
ymin = min(History.costDifferences);
ymax = max(History.costDifferences);
if ymin ~= ymax, ylim([ymin ymax]); end
xlim([1 nIter])
box on
clear h   
set(gcf,'Color',[1 1 1])
set(gcf,'PaperPositionMode','auto','InvertHardcopy','off')
set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.8])
figPath = fullfile(figDir,'Cost Difference History');
print(strcat(figPath,'.png'),'-dpng','-r250');
savefig(strcat(figPath,'.fig'))
close(gcf)

% Cost per Solution
for n = 1:nVariables
    figure
    hold on
    h(1) = scatter(History.solutions(n,iReject),History.costs(iReject),...
        20,[0.85 0.85 0.85],'filled');
    h(2) = scatter(History.solutions(n,iSuccess),History.costs(iSuccess),...
        20,[0.5 0.5 0.5],'filled');
    h(3) = scatter(History.solutions(n,iAccept),History.costs(iAccept),...
        20,[0 0 0],'filled');    
    title(sprintf('\\bfCost vs Variable %d \\rm(Opt. Cost = %0.3g)',...
        n,History.bestCost))
    xlabel(sprintf('Value of Variable x(%d)',n),'FontSize',10)
    ylabel('Cost C','FontSize',10)
    pbaspect([1 1 1])
    legend(h,{'rejected','accepted ("bad")','accepted ("good")'},...
        'Location','NorthEast','FontSize',10)
    xmin = min(History.solutions(n,:));
    xmax = max(History.solutions(n,:));
    ymin = min(History.costs);
    ymax = max(History.costs);
%     xmin = History.bestSolution(n) - std(History.solutions(n,:));
%     xmax = History.bestSolution(n) + std(History.solutions(n,:));
%     ymin = History.bestCost - 1e-3*std(History.costs);
%     ymax = History.bestCost + 5e-3*std(History.costs);
    if xmin ~= xmax, xlim([xmin xmax]); end
    if ymin ~= ymax, ylim([ymin ymax]); end
    box on
    clear h
    set(gcf,'Color',[1 1 1])
    set(gcf,'PaperPositionMode','auto','InvertHardcopy','off')
    set(gcf,'units','normalized','outerposition',[0.3 0.1 0.4 0.8])
    figPath = fullfile(figDir,sprintf('Cost vs Solution (Var%d)',n));
    print(strcat(figPath,'.png'),'-dpng','-r250');
    savefig(strcat(figPath,'.fig'))
    close(gcf)
end