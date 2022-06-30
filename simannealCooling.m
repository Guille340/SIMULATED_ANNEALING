%  To = simannealCooling(OptimValues,Problem)
%
%  DESCRIPTION
%  Cooling schedule for the computation of the temperature at the current
%  temperature step. The temperature step is given in OPTIMVALUES.TEMPITER.
%
%  There are four options available for the cooling schedule, given by the
%  'CoolingFcn' property in the OPTIONS structure: (1) 'coolingexp', for 
%  exponential cooling of the form TO = ALPHA^J * T0; (2) 'coolinglin', for 
%  linear cooling of the form TO = (1 - J/NTEMPSTEPS) * T0; (3) 'coolinginv'
%  for cooling with the inverse function, as in TO = T0/J; and (4) 'coolinglog' 
%  for optimal logarithmic cooling of the form TO = T0/(1 + LOG(J)). In these 
%  expressions, J is the temperature step, T0 is the initial temperature, ALPHA 
%  is the cooling factor (0 to 1), and TO is the temperature for the selected 
%  step.
%
%  SIMANNEALCOOLING receives the structure of optimisation parameters for the
%  current iteration (OPTIMVALUES) and the PROBLEM structure containing the
%  input data for SIMANNEAL.
%
%  A custom cooling schedule can be provided as a function handle through
%  the option 'CoolingFcn' in OPTIONS structure (see SIMANNEAL.m and 
%  SIMANNEALOPTIONS.m).
%
%  INPUT ARGUMENTS
%  - OptimValues: structure containing the optimisation parameters for the
%    current iteration. For details see INITIALISEOPTIMVALUES.m.
%  - Problem: structure containing the input settings used by SIMANNEAL.m.
%    # costFunction: objective function (FUN).
%    # initialSolution: initial multivariate solution (X0)
%    # lowerBoundaries: lower multivariate boundaries of solution space (LB).
%    # upperBoundaries: upper multivariate boundaries of solution space (UB).
%    # Options: options structure. For details, see SIMANNEALOPTIONS.m.
%
%  OUTPUT ARGUMENTS
%  - To: temperature for the specified step
%
%  FUNCTION CALL
%  1. To = simannealCooling(OptimValues,Problem)
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

function To = simannealCooling(OptimValues,Problem)

% Load Parameters
coolingFcn = Problem.Options.coolingFcn;
alpha = Problem.Options.coolingFactor;
T0 = Problem.Options.initialTemperature;
j = OptimValues.tempIter;
nTempSteps = OptimValues.nTempSteps;

% Process New Temperature
if isa(coolingFcn,'function_handle')
    To = coolingFcn(OptimValues,Problem);
else
    switch coolingFcn
        case 'coolingexp'
            To = alpha^j * T0;
        case 'coolinglin'
            To = (1 - j/nTempSteps) * T0;
        case 'coolinginv'
            To = T0/j;
        case 'coolinglog'
            To = T0/(1 + log(j));
    end
end

