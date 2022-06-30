%  OptimValues = initialiseOptimValues(Problem)
%
%  DESCRIPTION
%  Initialises the OPTIMVALUES structure. The structure is used as an internal
%  variable in SIMANNEAL.m to provide the various internal function with the 
%  values of the optimisation parameters for each iteration.
%
%  INITIALISEOPTIMVALUES takes the PROBLEM structure containing the input
%  settings for the problem.
%
%  INPUT ARGUMENTS
%  - PROBLEM: structure containing the inputs settings used by SIMANNEAL.m.
%    # costFunction: objective function (FUN).
%    # initialSolution: initial multivariate solution (X0)
%    # lowerBoundaries: lower multivariate boundaries of solution space (LB).
%    # upperBoundaries: upper multivariate boundaries of solution space (UB).
%    # Options: options structure. For details, see SIMANNEALOPTIONS.m.
%
%  OUTPUT ARGUMENTS
%  - OptimValues: structure containing the optimisation parameters for the 
%    current annealing iteration.
%    # solution: multivariate solution for current iteration.
%    # cost: cost of solution for current iteration.
%    # costDifference: difference in cost between the current and previous
%      colution.
%    # temperature: temperature step for the current iteration. Note that
%      the temperature typically remains constant for a number of iterations.
%      See property 'MaxIterPerTemp' in OPTIONS structure (SIMANNEALOPTIONS.m).
%    # globalIter: index of current iteration.
%    # tempIter: index of current temperature step.
%    # meanTemperature: mean temperature calculated over all temperature
%      iterations.
%    # nTempSteps: number of temperature steps. Calculated using the selected
%      cooling procedure, the initial temperature and the final temperature
%      (see 'CoolingFcn', 'InitialTemperature', 'minTemperature' properties
%      in OPTIONS structure).
%
%  FUNCTION CALL
%  1. OptimValues = initialiseOptimValues(Problem)
%
%  FUNCTION DEPENDENCIES
%  - None
%
%  TOOLBOX DEPENDENCIES
%  - MATLAB (Core)
%
%  See also SIMANNEAL

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  21 Jun 2022

function OptimValues = initialiseOptimValues(Problem)

% Load Parameters
x0 = Problem.initialSolution;
fun = Problem.costFunction;
T0 = Problem.Options.initialTemperature;
Tmin = Problem.Options.minTemperature;

% Initialise OptimValues
OptimValues.solution = x0;
OptimValues.cost = fun(x0);
OptimValues.costDifference = 0;
OptimValues.temperature = T0;
OptimValues.globalIter = 1;
OptimValues.tempIter = 1;
OptimValues.meanTemperature = NaN;
OptimValues.nTempSteps = NaN;

% Temporal OptimValues Structure
OptimValues_loop = OptimValues;

if ~isempty(T0)
    % Temperature Vector
    i = 1;
    Ti(i) = T0;

    % Calculate Number of Temperatures
    while OptimValues_loop.temperature > Tmin ...
            || OptimValues_loop.tempIter == 1e6
        % Calculate New Temperature
        Ti_new = simannealCooling(OptimValues_loop,Problem);

        % Update OptimValues
        OptimValues_loop.temperature = Ti_new;
        OptimValues_loop.tempIter = OptimValues_loop.tempIter + 1;

        % Temperature Vector
        i = i + 1;
        Ti(i) = Ti_new;
    end
    nTempSteps = OptimValues_loop.tempIter;

    % Update OptimValues
    OptimValues.meanTemperature = mean(Ti);
    OptimValues.nTempSteps = nTempSteps;
end
