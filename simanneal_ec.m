function [fun,x0,lb,ub,Options] = simanneal_ec(fun,x0,varargin_cell)

nFixArg = 2;
varargin = varargin_cell;
nVarargin = length(varargin);
nargin = nFixArg + nVarargin;

% Parse Input Arguments
switch nargin
    case 2
        lb = [];
        ub = [];
        Options = simannealOptions();
    case 3
        error('Wrong number of input arguments')
    case 4
        lb = varargin{1};
        ub = varargin{2};
        Options = simannealOptions();
    case 5
        lb = varargin{1};
        ub = varargin{2};
        Options = varargin{3};
end

% Error Control ('fun')
if ~isa(fun,'function_handle')
    error('FUN must be a function handle')
end

% Error Control ('x0')
if ~isnumeric(x0) || ~isvector(x0)
    error('X0 must be a numeric vector')
end
x0 = x0(:);

% Error Control ('lb')
if ~isempty(lb) && (~isnumeric(lb) || ~isvector(lb) || numel(lb) ~= numel(x0))
    error('LB must be a numeric vector the same size as X0')
end
lb = lb(:);

% Error Control ('ub')
if ~isempty(ub) && (~isnumeric(ub) || ~isvector(ub) || numel(ub) ~= numel(x0))
    error('UB must be a numeric vector the same size as X0')
end
ub = ub(:);

% Error Control ('Options')
if ~isOptionsStruct(Options)
    Options = initialiseOptionsStruct;
    warning(['OPTIONS is not a valid options structure. The default '...
        'values will be used (see INITIALISEOPTIONSSTRUCT for details '...
        'about the initalisation parameters)'])
end