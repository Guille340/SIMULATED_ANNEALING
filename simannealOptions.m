%  Options = simannealOptions(varargin)
%
%  DESCRIPTION
%  Creates a OPTIONS structure containing the settings for the simulated
%  annealing problem. OPTIONS is used as input for SIMANNEAL.m.
%
%  INPUT PROPERTIES
%  - 'AcceptanceFcn': function handle or character vector specifying the
%    the function to be used to accept or reject a newly generated test
%    solution. The function handle must be of the form @myfun(XI_NEW,CI_NEW,
%    OPTIMVALUES,PROBLEM). There is only one available character vector to
%    choose from, 'acceptancesa', which applies the Metropolis algorithm
%    (Metropolis et al., 1953).
%  - 'GeneratorFcn': function handle or character vector specifying the
%    the function to be used to perturb a test solution to generate a new
%    one. The function handle must be of the form @myfun(OPTIMVALUES,PROBLEM).
%    There is are two character vectors to choose from
%    # 'GeneratorVar': perturbation algorithm with the size of the perturbation 
%      step proportional to the temperature (default).
%    # GeneratorCst': perturbation algorithm with a constant perturbation step.
%  - 'CoolingFcn': function handle or character vector specifying the
%    the cooling function to be used to calculate the temperature for a given
%    temperature step. The function handle must be of the form @myfun(...
%    OPTIMVALUES,PROBLEM). There is are four character vectors to choose from
%    # 'coolingexp': exponential cooling of the form TO = ALPHA^J * T0 (default)
%    # 'coolinglin': linear cooling of the form TO = (1 - J/NTEMPSTEPS) * T0
%    # 'coolinginv': cooling with the inverse function, as in TO = T0/J
%    # 'coolinglog': for optimal logarithmic cooling of the form TO = 
%      T0/(1 + LOG(J))
%    In the above expressions, J is the temperature step, T0 is the initial 
%    temperature, ALPHA is the cooling factor (0 to 1), and TO is the 
%    temperature for the selected step.
%  - 'InitialTemperature': initial (maximum) temperature. This parameter is
%    critical for the correct annealing. In general, a value comparable to
%    the mean positive cost difference of a few randomly generated transitions 
%    is a good start. The value can also be calculate automatically by setting
%    'InitialTemperature' = [] (default).
%  - 'CoolingFactor': value between 0 and 1 used for the exponential cooling
%    schedule (TO = COOLINGFACTOR^J * T0). COOLINGFACTOR = 0.9 (default).
%  - 'MaxIterations': maximum total number of iterations. MAXITERATIONS = INF
%    (default).
%  - 'MaxIterPerTemp': maximum number of iterations within a temperature.
%     MAXITERPERTEMP = 300 (default).
%  - 'MaxAccept': maximum number of accepted transitions within a temperature.
%     MAXACCEPT = 20 (default).
%  - 'MaxConsReject': maximum number of consecutive rejected transitions 
%     within a temperature. MAXCONSREJECT = 1000 (default).
%  - 'MinTemperature': minimum temperature. 'MinTemperature', 'CoolingFcn' and
%    'InitialTemperature' determine how fast the annealing process occurs
%     (i.e., number of temperature steps). A larger ratio between the initial
%     (maximum) and minimum temperature will make the process slower (and in
%     many cases, more accurate). MINTEMPERATURE = 1e-8
%  - 'Display': character vector specifying the type of report to provide on
%     the command window.
%     # 'final': displays one final report (default).
%     # 'iter': displays one report per temperature step and a final report.
%     # 'off': disabled.
%  - 'Plot': character vector indicating if the results are to be plotted.
%     # 'on': plot results
%     # 'off': do not plot results.
%
%     NOTE: 'MaxIterations', 'MaxIterPerTemp', 'MaxAccept', 'MaxConsReject',
%     and 'MinTemperature' are stop criteria. If these numbers are exceeded
%     at any point during the execution, the algorithm stops and returns the 
%     best solution, best cost and History up to this point.
%
%  OUTPUT ARGUMENTS
%  - OPTIONS: structure containing all the available configuration fields.
%    Those not specified as input are included with the default values. The
%    fields have the following names: acceptanceFcn, generatorFcn, coolingFcn,
%    coolingFcn, initialTemperature, coolingFactor, maxIterations, 
%    maxIterPerTemp, maxAccept, maxConsReject, minTemperature, display, and
%    plot. Their meaning can be deduced from the description of the input
%    properties.
%
%  FUNCTION CALL
%  1. Options = simannealOptions(PROPERTY,VALUE)
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

function Options = simannealOptions(varargin)

% Verify number of Input Arguments
nFixArg = 0;
nVarargin = length(varargin);
nargin = nFixArg + nVarargin;
if rem(nargin - nFixArg,2)
    error('Variable input arguments must come in pairs (PROPERTY,VALUE)')
end

% Extract and Verify Input Properties
validProperties = {'acceptancefcn','generatorfcn','coolingfcn',...
    'initialtemperature','coolingfactor','maxiterations',...
    'maxiterpertemp','maxaccept','maxconsreject','mintemperature',...
    'display','plot'};
properties = lower(varargin(1:2:end));
if any(~ismember(properties,validProperties))
    error('One or more input properties are not recognised')
end

% Default Input Values
acceptanceFcn = 'acceptancesa';
generatorFcn = 'generatorvar';
coolingFcn = 'coolingexp';
T0 = [];
alpha = 0.9;
maxIter = Inf;
maxIterPerTemp = 300;
maxAccept = 20;
maxConsReject = 1000;
Tmin = 1e-8;
displayOpt = 'final';
plotOpt = 'off';

% Extract and Verify Input Values
values = varargin(2:2:end);
nPairs = (nargin - nFixArg)/2; % number of (PROPERTY,VALUE) pairs
for m = 1:nPairs
    property = properties{m};
    switch property % populate with more properties if needed
        case 'acceptancefcn'
            acceptanceFcn = values{m};
            if ~isa(acceptanceFcn,'function_handle') ...
                    && (~ischar(acceptanceFcn) || ~strcmpi(acceptanceFcn,...
                    'acceptancesa'))
                acceptanceFcn = 'acceptancesa';
                warning(['Non-supported value for PROPERTY = '...
                    '''AcceptanceFcn''. The input must be a valid function '...
                    'handle or supported character string. ACCEPTANCEFCN = '...
                    '''%s'' will be used'],acceptanceFcn)  
            end
        case 'generatorfcn'
            generatorFcn = values{m};
            if ~isa(generatorFcn,'function_handle') ...
                    && (~ischar(generatorFcn) ...
                    || ~any(strcmpi(generatorFcn,...
                    {'generatorvar','generatorcst'})))
                generatorFcn = 'generatorvar';
                warning(['Non-supported value for PROPERTY = '...
                    '''GeneratorFcn''. The input must be a valid function '...
                    'handle or supported character string. GENERATORFCN = '...
                    '''%s'' will be used'],generatorFcn)  
            end
        case 'coolingfcn'
            coolingFcn = values{m};
            if ~isa(coolingFcn,'function_handle') ...
                    && (~ischar(coolingFcn) ...
                    || ~any(strcmpi(coolingFcn,...
                    {'coolingexp','coolinglin','coolinginv','coolinglog'})))
                coolingFcn = 'coolingexp';
                warning(['Non-supported value for PROPERTY = '...
                    '''CoolingFcn''. The input must be a valid function '...
                    'handle or supported character string. COOLINGFCN = '...
                    '''%s'' will be used'],coolingFcn)  
            end
        case 'initialtemperature'
            T0 = values{m};
            if ~isempty(T0) && (~isnumeric(T0) || ~isscalar(T0) ...
                    || T0 < 0)
                T0 = [];
                warning(['Non-supported value for PROPERTY = '...
                    '''InitialTemperature''. The input must be a positive '...
                    'scalar. The initial temperature will be '...
                    'calculated automatically (INITIALTEMPERATURE = [])'])
            end
        case 'coolingfactor'
            alpha = values{m};
            if ~isnumeric(alpha) || ~isscalar(alpha) || alpha < 0 || alpha > 1
                alpha = 0.9;
                warning(['Non-supported value for PROPERTY = '...
                    '''CoolingFactor''. The input must be a positive '...
                    'scalar between 0 and 1. COOLINGFACTOR = %0.1f will '...
                    'be used'],alpha)
            end
        case 'maxiterations'
            maxIter = values{m};
            if ~isnumeric(maxIter) || ~isscalar(maxIter) || rem(maxIter,1) ...
                    || maxIter < 0
                maxIter = Inf;
                warning(['Non-supported value for PROPERTY = '...
                    '''MaxIterations''. The input must be a positive '...
                    'integer. MAXITERATIONS = %0.1f will be used'],maxIter)
            end
        case 'maxiterpertemp'
            maxIterPerTemp = values{m};
            if ~isnumeric(maxIterPerTemp) || ~isscalar(maxIterPerTemp) ...
                    || rem(maxIterPerTemp,1) || maxIterPerTemp < 0
                maxIterPerTemp = 300;
                warning(['Non-supported value for PROPERTY = '...
                    '''MaxIterPerTemp''. The input must be a positive '...
                    'integer. MAXITERATIONSPERTEMP = %0.1f will be used'],...
                    maxIterPerTemp)
            end
        case 'maxaccept'
            maxAccept = values{m};
            if ~isnumeric(maxAccept) || ~isscalar(maxAccept) ...
                    || rem(maxAccept,1) || maxAccept < 0
                maxAccept = 20;
                warning(['Non-supported value for PROPERTY = '...
                    '''MaxAccept''. The input must be a positive '...
                    'integer. MAXACCEPT = %0.1f will be used'],maxAccept)
            end
        case 'maxconsreject'
            maxConsReject = values{m};
            if ~isnumeric(maxConsReject) || ~isscalar(maxConsReject) ...
                    || rem(maxConsReject,1) || maxConsReject < 0
                maxConsReject = 1000;
                warning(['Non-supported value for PROPERTY = '...
                    '''MaxConsReject''. The input must be a positive '...
                    'integer. MAXCONSREJECT = %0.1f will be used'],...
                    maxConsReject)
            end
        case 'mintemperature'
            Tmin = values{m};
            if ~isnumeric(Tmin) || ~isscalar(Tmin) || Tmin < 0
                Tmin = 1e-8;
                warning(['Non-supported value for PROPERTY = '...
                    '''MinTemperature''. The input must be a positive '...
                    'number. MINTEMPERATURE = %0.1f will be used'],...
                    Tmin)
            end
        case 'display'
            displayOpt = values{m};
            if ~ischar(displayOpt) || ~any(strcmpi(displayOpt,...
                    {'off','iter','final'}))
                displayOpt = 'final';
                warning(['Non-supported value for PROPERTY = '...
                    '''Display''. The input must be a supported character '...
                    'string. DISPLAY = ''%s'' will be used'],displayOpt)  
            end
        case 'plot'
            plotOpt = values{m};
            if ~ischar(plotOpt) || ~any(strcmpi(plotOpt,{'off','on'}))
                plotOpt = 'off';
                warning(['Non-supported value for PROPERTY = '...
                    '''Plot''. The input must be a supported character '...
                    'string. PLOT = ''%s'' will be used'],plotOpt)  
            end
    end
end

% Build Output Structure
Options.acceptanceFcn = acceptanceFcn;
Options.generatorFcn = generatorFcn;
Options.coolingFcn = coolingFcn;
Options.initialTemperature = T0;
Options.coolingFactor = alpha;
Options.maxIterations = maxIter;
Options.maxIterPerTemp = maxIterPerTemp;
Options.maxAccept = maxAccept; 
Options.maxConsReject = maxConsReject;
Options.minTemperature = Tmin;
Options.display = displayOpt;
Options.plot = plotOpt;
