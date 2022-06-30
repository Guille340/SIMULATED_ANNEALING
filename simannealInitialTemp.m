%  T0 = simannealInitialTemp(Problem,varargin)
%
%  DESCRIPTION
%  Calculates the initial temperature from a few randomly generated transitions
%  (T -> Inf). The calculation METHOD and number of iterations (NITER) can be 
%  specified.
%
%  INPUT ARGUMENTS
%  - Problem: structure containing the input settings used by SIMANNEAL.m.
%    # costFunction: objective function (FUN).
%    # initialSolution: initial multivariate solution (X0)
%    # lowerBoundaries: lower multivariate boundaries of solution space (LB).
%    # upperBoundaries: upper multivariate boundaries of solution space (UB).
%    # Options: options structure. For details, see SIMANNEALOPTIONS.m.
%  - method (varargin{1}): calculation method
%    # 'johnson': T0 = -mean(DC)/log(Pa0), where DC is the cost difference
%      from a set of randomly generated transitions and PA0 is an initial 
%      probability of acceptance estimated as PA0 = 0.8. This is the default
%      option. [Johnson et al., 1989; Xinchao, 2011]
%    # 'aarts': T0 = K*std(C)^2, where STD(C) is the standard deviation of the
%       costs from a set of randomly generated transitions and K is a constant
%       between 5 and 10 (K = 8 is assummed) [Aarts et al, 1997]
%    # 'kirkpatrick': T0 = MAX(C) - MAX(C), where MAX(C) and MIN(C) are the 
%       maximum and minimum cost of a set of randomly generated transitions.
%  - nIter: number of randomly generated transitions. NITER = 100 (default)
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

function T0 = simannealInitialTemp(Problem,varargin)

narginchk(1,3)

if nargin == 1
    method = 'johnson';
    nIter = 100;
elseif nargin == 2
    method = varargin{1};
    nIter = 100;
else
    method = varargin{1};
    nIter = varargin{2};
end

% LOAD PARAMETERS
fun = Problem.costFunction;
       
% CALCULATE COST ON A NUMBER OF RANDOM SOLUTIONS (T0 -> Inf)
% Initialise OptimValues
OptimValues = initialiseOptimValues(Problem);

% Initialise T and T0 to Identical Finite Numbers
OptimValues.temperature = 1;
Problem.Options.initialTemperature = 1; 

% Initialise Cost Vector
C = OptimValues.cost;

% Repeat for All Iterations
for i = 2:nIter % 100 iterations the cost seems to stabilise for most cases
        
    % Compute New Parameters
    xi_new = simannealGenerator(OptimValues,Problem);
    Ci_new = fun(xi_new);
    
    % Update Cost Vector
    C(i) = Ci_new;
    
    % Update OptimValues
    OptimValues.solution = xi_new;
    OptimValues.cost = Ci_new;
    OptimValues.globalIter = i;
    OptimValues.tempIter = OptimValues.tempIter + 1;
end

switch lower(method)
    case 'johnson' % [Johnson et al., 1989; Xinchao, 2011]
        Pa0 = 0.8; % initial probability of acceptance
        DC = diff(C);
        DC = DC(DC >= 0); % use only positive cost differences ("bad" jumps)
        T0 = -mean(DC)/log(Pa0);
        
    case 'aarts' % [Aarts et al, 1997]
        K = 8; % constant (typ. between 5-10)
        T0 = K*std(C)^2;
        
    case 'kirkpatrick' % [Kirpatrick, 1983]
        T0 = max(C) - min(C);
end

