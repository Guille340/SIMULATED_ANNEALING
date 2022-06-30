%  [accept,P,r] = simannealAcceptance(xi_new,Ci_new,OptimValues,Problem)
%
%  DESCRIPTION
%  Acceptance criteria for the selection or rejection of an input test
%  solution XI_NEW with cost CI_NEW. The function receives the structure 
%  OPTIMVALUES containing the optimisation parameters for the current iteration
%  and the PROBLEM structure containing the input settings for the annealing
%  process. For details about the content of the OPTIMVALUES structure, see 
%  function INITIALISEOPTIMVALUES.m.
% 
%  SIMANNEALACCEPTANCE.m returns the acceptance flag, probability of acceptance
%  and the uniform random number used to evaluate that probability. The 
%  function applies the Metropolis algorithm to determine the acceptance or
%  rejection of test solution XI_NEW (Metropolis et al., 1953). The temperature
%  and the previous solution and its cost are included in OPTIMVALUES.
%
%  A custom acceptance schedule can be provided as a function handle through
%  the option 'AcceptanceFcn' in OPTIONS structure (see SIMANNEAL.m and 
%  SIMANNEALOPTIONS.m).
%
%  INPUT ARGUMENTS
%  - XI_NEW: new multivariate test solution.
%  - CI_NEW: cost of new multivariate test solution.
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
%  - accept: acceptance flag. There are three possible values.
%    # 1: solution from "good" jump (COST_NEW - COST_OLD < 0) is always 
%      accepted .
%    # -1: solution from "bad" jump (COST_NEW - COST_OLD > 0) is accepted 
%      if the probability of acceptance is higher than a random number 
%      between 0 and 1 (P > R)
%    # 0: solution from "bad" jump (COST_NEW - COST_OLD > 0) is rejected if 
%      the probability of acceptance is lower than a random number 
%      between 0 and 1 (P < R).
%  - P: probability of acceptance.
%  - r: uniform random number between 0 and 1 generated to determine whether
%    a solution XI_NEW resulting in a "bad" jump is to be accepted (P > r)
%    or rejected (P <= r).
%
%  FUNCTION CALL
%  1. accept = simannealAcceptance(xi_new,OptimValues,Problem)
%  2. [accept,P,r] = simannealAcceptance(...)
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

function [accept,P,r] = simannealAcceptance(xi_new,Ci_new,OptimValues,Problem)

% Load Parameters
acceptanceFcn = Problem.Options.acceptanceFcn;
Ci = OptimValues.cost;
Ti = OptimValues.temperature;

% Run Acceptance Function
if isa(acceptanceFcn,'function_handle')
    [accept,P,r] = acceptanceFcn(xi_new,Ci_new,OptimValues,Problem);
else % acceptanceFcn = 'acceptancesa'
    % Compute Cost
    DC = Ci_new - Ci;

    % Apply Acceptance Criteria
    if DC <= 0 % direct acceptance of solution ("good jump")
        accept = 1; % accepted "new improved solution"
        P = 1;
        r = NaN;
    else % Metropolis algorithm ("bad jump")
        P = exp(-DC/Ti); % probability of acceptance of new solution
        r = rand;
        if P > r
            accept = -1; % accepted "new less-optimal solution"
        else
            accept = 0; % rejected "new less-optimal solution"
        end
    end
end
