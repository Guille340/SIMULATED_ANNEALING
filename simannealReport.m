%  simannealReport(History,reportType)
%
%  DESCRIPTION
%  Displays a report of the results on the command window. The details of the
%  report are determined by the REPORTTYPE option. The function uses the
%  HISTORY structure returned from SIMANNEAL.m to gather all the necessary
%  information.
%
%  INPUT ARGUMENTS
%  - History: structure containing the optimisation parameters for every 
%    iteration happening in the annealing process. Generated with function
%    SIMANNEAL.m.
%  - reportType: this is one of the character vectors available for the 
%    'Display' option in the OPTIONS structure.
%    # 'final': displays one final report (default).
%    # 'iter': displays one report per temperature step and a final report.
%    # 'off': disabled.
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

function simannealReport(History,reportType)

switch reportType
    case 'final'
        % Calculate Parameters
        x_best = History.bestSolution;
        C_best = History.bestCost;
        C_mean = mean(History.costs);
        C_std = std(History.costs);
        T0 = History.temperatures(1);
        T_last = History.temperatures(end);
        nTempSteps = History.nTempSteps;
        nIter = length(History.costs);
        nAccept = sum(History.accepts == 1);
        nSuccess = sum(History.accepts == -1);
        nReject = sum(History.accepts == 0);
        isReject = diff([0 History.accepts == 0 0]);
        nConsReject = max(find(isReject == -1) - find(isReject == 1));
        
        % Report
        fprintf('\n\n [Final Report]')
        fprintf('\n Solution = ')
        fprintf('%0.3g ',x_best) % solution
        fprintf('\n Cost = %0.3g',C_best) % cost
        fprintf('\n Cost Mean = %0.3g',C_mean) % cost average
        fprintf('\n Cost StdDev = %0.3g',C_std) % cost standard deviation
        fprintf('\n Initial Temperature = %0.3g',T0) % initial temperature
        fprintf('\n Final Temperature = %0.3g',T_last) % initial temperature
        fprintf('\n No. Temperature Steps = %d',nTempSteps) % number of temperature steps
        fprintf('\n No. Iterations = %d',nIter) % number of iterations
        fprintf('\n No. Accepted ("good jump") = %d',nAccept) % number of accepted 'good' solutions
        fprintf('\n No. Accepted ("bad jump") = %d',nSuccess) % number of accepted 'bad' solutions
        fprintf('\n No. Rejected = %d',nReject) % number of rejected 'bad' solutions
        fprintf('\n No. Rejected (Consec.) = %d',nConsReject) % number of rejected 'bad' solutions
        fprintf('\n Time = %s\n\n',datestr(clock)) % total elapsed time (since beginning of calculations) 
        
    case 'iter'
        
        % Slect Elements from Last Temperature
        tempStep = length(unique(History.temperatures)); % no. last temperature step
        ind = find(History.temperatures == History.temperatures(end));
        x = History.solutions(:,ind);
        C = History.costs(ind);
        T = History.temperatures(ind);
        accepts = History.accepts(ind);
        
        % Paramters from Last Temperature
        [~,iBest] = min(C);
        C_best = C(iBest);
        x_best = x(:,iBest);
        T_best = T(end);
        C_mean = mean(C);
        C_std = std(C);
        nIter = length(ind);
        nAccept = sum(accepts == 1);
        nSuccess = sum(accepts == -1);
        nReject = sum(accepts == 0);
        isReject = diff([0 accepts == 0 0]);
        nConsReject = max(find(isReject == -1) - find(isReject == 1));
        if isempty(nConsReject), nConsReject = 0; end
        
        % Report
        fprintf('\n\n [Partial Report (%d/%d)]',tempStep,History.nTempSteps)
        fprintf('\n Temperature = %0.3g',T_best) % current temperature
        fprintf('\n Solution = ')
        fprintf('%0.2g ',x_best) % last solution within temperature
        fprintf('\n Cost = %0.3g',C_best) % last cost within temperature
        fprintf('\n Cost Mean = %0.3g',C_mean) % cost average
        fprintf('\n Cost StdDev = %0.3g',C_std) % cost standard deviation
        fprintf('\n No. Iterations = %d',nIter) % number of iterations within temperature
        fprintf('\n Accepted ("good jump") = %d',nAccept) % number of accepted 'good' solutions
        fprintf('\n Accepted ("bad jump") = %d',nSuccess) % number of accepted 'bad' solutions
        fprintf('\n No. Rejected = %d',nReject) % number of rejected 'bad' solutions
        fprintf('\n No. Rejected (Consec.) = %d',nConsReject) % number of rejected 'bad' solutions
        fprintf('\n Time = %s\n\n',datestr(clock)) % total elapsed time (since beginning of calculations) 

    case 'off'
           
end
