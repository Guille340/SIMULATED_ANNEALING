%  xo = simannealGenerator(OptimValues,Problem)
%
%  DESCRIPTION
%  Perturbation algorithm for the calculation of a new test solution through
%  a small perturbation applied to the current solution. The current solution
%  is included in the structure OPTIMVALUES.
%
%  There are two options available for the pertubation algorithm, given by the
%  'GeneratorFcn' property in the OPTIONS structure: (1) 'GeneratorVar', for
%  a pertubation step whose size depends on the temperature, (2) 'GeneratorCst',
%  for a constant perturbation step.
%
%  A custom perturbation algorithm can be provided as a function handle through
%  the option 'GeneratorFcn' in OPTIONS structure (see SIMANNEAL.m and 
%  SIMANNEALOPTIONS.m).
%
%  INPUT ARGUMENTS
%  - OPTIMVALUES: structure containing the optimisation parameters for the
%    current iteration. For details see INITIALISEOPTIMVALUES.m.
%  - PROBLEM: structure containing the input settings used by SIMANNEAL.m.
%    # costFunction: objective function (FUN).
%    # initialSolution: initial multivariate solution (X0)
%    # lowerBoundaries: lower multivariate boundaries of solution space (LB).
%    # upperBoundaries: upper multivariate boundaries of solution space (UB).
%    # Options: options structure. For details, see SIMANNEALOPTIONS.m.
%
%  OUTPUT ARGUMENTS
%  - xo: new test solution obtained through perturbation.
%
%  FUNCTION CALL
%  1. xo = simannealGenerator(OptimValues,Problem)
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

function xo = simannealGenerator(OptimValues,Problem)

% Load Parameters
lb = Problem.lowerBoundaries;
ub = Problem.upperBoundaries;
generatorFcn = Problem.Options.generatorFcn;
T0 = Problem.Options.initialTemperature;
xi = OptimValues.solution;
Ti = OptimValues.temperature;
T_mean = OptimValues.meanTemperature;

if isa(generatorFcn,'function_handle')
    xo = generatorFcn(OptimValues,Problem);
else
    switch generatorFcn
        case 'generatorvar'
            T = Ti;
        case 'generatorcst'
            T = T_mean;
    end
    
    % OPTION 1: PERTURB ONE VARIABLE
    % Logical Vector of Perturbed Variable
    iDyn = find(lb ~= ub); % TRUE for non-constant variables
    iPerturb = iDyn(randperm(length(iDyn),1)); % perturb only non-constant variables
    isVar = false(length(xi),1);
    isVar(iPerturb) = true; % variable to be perturbed is set as TRUE

    % Step Size
    emptyBounds = isempty(lb) || any(isinf(lb)) || isempty(ub) || any(isinf(ub));
    normFactor = 1;
    if ~emptyBounds
        normFactor = 1e0 * abs(ub(isVar) - lb(isVar))/T0; % normalise by boundary size
    end
    stepSize = T * randn * normFactor; % max step amplitude ~T*normFactor

    % Compute Perturbed Solution
    xo = xi + isVar*stepSize;
    if ~emptyBounds
        endloop = (xo(isVar) >= lb(isVar)) && (xo(isVar) <= ub(isVar));
        while ~endloop % recalculate solution if exceeds solution space
            stepSize = T * randn * normFactor;
            xo = xi + isVar*stepSize;
            endloop = (xo(isVar) >= lb(isVar)) && (xo(isVar) <= ub(isVar));
        end
    end

    % OPTION 2: PERTURB ALL VARIABLES
    % % General
    % nVar = length(xi); % number of variables in input test solution xi
    % 
    % % Step Size
    % emptyBounds = isempty(lb) || isinf(lb) || isempty(ub) || isempty(ub);
    % if emptyBounds
    %     stepSize = T * (2*round(rand(nVar,1)) - 1); % +/- temperature
    % else
    %     stepSize = T/T0 * 0.2*abs(ub - lb) .* (2*round(rand(nVar,1)) - 1);
    % end
    % 
    % % Compute Perturbed Solution
    % xo = xi + stepSize;
    % if ~emptyBounds
    %     endloop = all(xo >= lb) && all(xo <= ub);
    %     while ~endloop % recalculate solution if exceeds solution space
    %         xo = xi + stepSize;
    %         endloop = all(xo >= lb) && all(xo <= ub);
    %     end
    % end
end
