%  [x,C,History] = simanneal(fun,x0,varargin)
%
%  DESCRIPTION
%  Finds the minimum value C and associated multivariate solution X of an 
%  input function FUN using the simulated annealing algorithm. The user
%  must provide a provisional solution X (start point). The lower and upper
%  boundaries of the problem (LB,UB) may also be provided. Specific settings
%  for the algorithm are given through an input OPTIONS structure generated
%  with function SIMANNEALOPTIONS.m. SIMANNEAL also returns the structure
%  HISTORY containing the evolution of the optimisation parameters, iteration-
%  by-iteration.
%
%  INPUT ARGUMENTS
%  - FUN: objective or cost function for which SIMANNEAL will attempt to
%    find the minimum C and associated optimal solution X.
%  - X0: initial multivariate solution. X0 must be within the limits of the
%    solution space (LB,UB), if provided.
%  - LB (varargin{1}): lower boundaries of the solution space. This is a 
%    vector with the same number of elements as X0.
%  - UB (varargin{2}): upper boundaries of the solution space. This is a
%    vector with the same number of elements as X0.
%  - OPTIONS: structure containing the settings for the simulated annealing 
%    algorithm. For details about its content and how to generate it, see
%    function SIMANNEALOPTIONS.m.
%
%  OUTPUT ARGUMENTS
%  - x: best solution found during simulated annealing.
%  - C: minimum cost associated with the best solution X. This is the minimum
%    value obtained for input function FUN.
%  - History: structure containing the optimisation parameters used during
%    simulated annealing process.
%    # bestSolution: best multivariate solution found during annealing.
%    # bestCost: minimum value of input cost function found during annealing.
%    # solutions: matrix of solutions obtained during annealing. One solution
%      per iteration. The solutions are organised in columns (matrix of 
%      dimensions [NSOLUTIONS,NITERATIONS]).
%    # temperatures: vector of temperatures used in the annealing process.
%      One temperature per iteration. Note that several iterations will take 
%      place before the temperature is reduced. For details about the cooling 
%      schedule see function SIMANNEALCOOLING.m.
%    # accepts: vector of acceptance flags. One value per iteration. For 
%      details see SIMANNEALACCEPTANCE.m . There are three possible values.
%      ¬  1: solution from "good" jump (COST_NEW - COST_OLD < 0) is always 
%        accepted .
%      ¬ -1: solution from "bad" jump (COST_NEW - COST_OLD > 0) is accepted 
%        if the probability of acceptance is higher than a random number 
%        between 0 and 1 (P > R)
%      ¬  0: solution from "bad" jump (COST_NEW - COST_OLD > 0) is rejected 
%        if the probability of acceptance is lower than a random number 
%        between 0 and 1 (P < R).
%    # probAccept: probability of acceptance of a "bad" jump (COST_NEW - 
%      COST_OLD > 0). One probability per iteration. 
%    # nTempSteps: number of temperature steps. Determined by the cooling
%      schedule and the minimum temperature. The minimum temperature is
%      determined with function SIMANNEALOPTIONS.m. For details about the
%      cooling schedule see function SIMANNEALCOOLING.m.
%
%  FUNCTION CALL
%  1. x = simanneal(fun,x0)
%  2. x = simanneal(fun,x0,lb,ub)
%  3. x = simanneal(fun,x0,lb,ub,Options)
%  4. [x,C] = simanneal(...)
%  5. [x,C,History] = simanneal(...)
%
%  FUNCTION DEPENDENCIES
%  - simanneal_ec
%  - simannealInitialTemp
%  - initialiseOptimValues
%  - simannealGenerator
%  - simannealCooling
%  - simannealAcceptance
%  - simannealReport
%  - simannealPlot
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  See also SIMANNEAL_EC

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Jun 2022

function [x,C,History] = simanneal(fun,x0,varargin)

narginchk(2,5)

% Error Control
[fun,x0,lb,ub,Options] = simanneal_ec(fun,x0,varargin);

% Build 'Problem' Structure
Problem.costFunction = fun;
Problem.initialSolution = x0;
Problem.lowerBoundaries = lb;
Problem.upperBoundaries = ub;
Problem.Options = Options;

% Parameters from 'Options' Structure
T0 = Options.initialTemperature; % initial temperature
maxIter = Options.maxIterations; % maximum number of iterations within temperature (stop criterium for temperature loop)
maxIterPerTemp = Options.maxIterPerTemp;
maxAccept = Options.maxAccept; % maximum number of accepted 'good' solutions within temperature (stop criterium for temperature loop)
maxConsReject = Options.maxConsReject; % maximum number of consecutive rejections (stop criterium for general loop)
Tmin = Options.minTemperature;
displayOpt = Options.display;
plotOpt = Options.plot;

% Calculate Initial Temperature (if not specified)
if isempty(T0)
    T0 = simannealInitialTemp(Problem);
    Problem.Options.initialTemperature = T0;
end

% Initialise OptimValues
OptimValues = initialiseOptimValues(Problem);

% Initialise History  
History.bestSolution = OptimValues.solution;
History.bestCost = OptimValues.cost;
History.solutions = OptimValues.solution; % solution (design vector)
History.costs = OptimValues.cost; % cost
History.temperatures = OptimValues.temperature; % temperature
History.accepts = 1;
History.probAccept = 1;
History.nTempSteps = OptimValues.nTempSteps;

% Simulated Annealing Algorithm (Initialise)
i = 1; % general iteration index
j = 1; % 1-temperature iteration index
cnt_temp = 1; % temperature step counter (T0 is considered the first step)
success = 0;
consrej = 0;
finished = 0;

% Simulated Annealing Algorithm (Process)
while ~finished

    % New Parameters
    xi_new = simannealGenerator(OptimValues,Problem);
    Ci_new = fun(xi_new);
    DCi_new = Ci_new - OptimValues.cost;
    Ti_new = OptimValues.temperature;

    % New Temperature
    if j >= maxIterPerTemp || success >= maxAccept
        % Report History (per temperature)
        if strcmp(displayOpt,'iter')
            simannealReport(History,'iter')     
        end
        
        % Calculate New Temperature
        Ti_new = simannealCooling(OptimValues,Problem); % reduce temperature one step
        cnt_temp = cnt_temp + 1; % increase temperature step counter by 1
        j = 1; % reset temperature loop counter
        success = 0; % reset good-perturbation counter
    end

    % Apply Acceptance Criteria
    [accept,P,r] = simannealAcceptance(xi_new,Ci_new,OptimValues,Problem);
    switch accept
        case 1 % accepted "new improved solution"
            success = success + 1;
            consrej = 0; % initialise consecutive rejections variable
        case -1 % accepted "new less-optimal solution"
            consrej = 0; % initialise consecutive rejections variable
        case 0 % rejected "new less-optimal solution"
            xi_new = OptimValues.solution; % keep old solution
            Ci_new = OptimValues.cost;
            consrej = consrej + 1; % increase consecutive rejections by 1 
    end

    % Stop Criteria (Temperature)
    if Ti_new < Tmin
        finished = true;
    end

    % Stop Criteria (Consecutive Rejections)
    if consrej > maxConsReject
        finished = true;
    end

    % Stop Criteria (Total Number of Iterations)
    if i == maxIter
        finished = true;
    end
    
    % Update Counters
    i = i + 1;
    j = j + 1;

    % Update OptimValues
    OptimValues.solution = xi_new;
    OptimValues.cost = Ci_new;
    OptimValues.costDifference = DCi_new;
    OptimValues.temperature = Ti_new;
    OptimValues.globalIter = i;
    OptimValues.tempIter = cnt_temp;

    % Update History
    History.solutions(:,i) = OptimValues.solution;
    History.costs(i) = OptimValues.cost;
    History.costDifferences(i) = OptimValues.costDifference;
    History.temperatures(i) = OptimValues.temperature;
    History.accepts(i) = accept;
    History.probAccept(i) = P;
    History.randomNum(i) = r;
    
    % Update Best Parameters
    [~,iBest] = min(History.costs);
    History.bestSolution = History.solutions(:,iBest);
    History.bestCost = History.costs(:,iBest);
end

% Output Variables
x = History.bestSolution;
C = History.bestCost;

% Report History (summary)
if any(strcmp(displayOpt,{'iter','final'}))
    simannealReport(History,'final')               
end

% Plot
if strcmpi(plotOpt,'on')
    simannealPlot(History)
end
