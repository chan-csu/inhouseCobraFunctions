function printModel(model, varargin)
% Print model constraints or variable bounds or everything, with the
% associated shadow price and reduced cost if a solution is supplied.
% Support problem structure for Gurobi, Cobra, or Cplex construct.
%
% USAGE:
%    printModel(model, itemsToPrint, typeToPrint, solToPrint, colNames, rowNames, fileToPrint, modelFormat)
%    printModel(model, 'parameter', value, ...)
%
% INPUT:
%    model:     Gurobi or Cobra optimization model structure that can be solved by
%               gurobi or solveCobraLP/MILP/QP, or Cplex Matlab construct.
%               Cobra model is also acceptable.
%
% OPTIONAL INPUTS:
%    (also support parameter-value inputs)
%    itemsToPrint:  index vector or cell array of strings for the constraints/variables to be printed.
%                   'all' for all items in `typeToPrint` (default 'all').
%    typeToPrint:   - 'sum' or 's' prints summary (default if `itemsToPrint' is empty or not supplied)
%                   - 'all' or 'a' prints the entire model (default if `itemsToPrint' is 'all')
%                   - 'obj' or 'o' prints the objective function
%                   - 'con' or 'c' prints constraints (default if `itemsToPrint` is an index vector)
%                   - 'var' or 'v' prints variable types and bounds
%                     Special constraints supported currently:
%                   - 'ind' print indicator constraints (Cplex only)
%                   - 'sos' print SOS constraints (Gurobi and Cplex only)
%    solToPrint:    a corresponding solution structure for the model. When supplied,
%                   will print also the variable values, reduced cost and shadow price 
%                   associated with the requested variables and constraints
%    colNames:      cell array of names for the variables. Used if not detected in model. Default {'x1', ..., 'xN'}
%    rowNames:      cell array of names for the constraints. Used if not detected in model. Default {'r1', ..., 'rN'}
%    fileToPrint:   name of the textfile for saving the print results. 1 or [] for printing on the command window (default)
%    modelFormat:   'gurobi', 'cplex' or 'cobra' to state which format the input problem is in (default automatically detected)

%% argument handling modified from `functionTemplate`

%%%%% add here other special items to print
otherPrintTypes = {'sos', 'ind'};
otherPrintTypesNames = struct('sos', 'SOS constraints', 'ind', 'indicator constraints');

optArgin = {     'itemsToPrint', 'typeToPrint', 'solToPrint', 'colNames', 'rowNames', 'fileToPrint', 'modelFormat'};
defaultValues = {NaN,             NaN,            [],           [],         [],         [],            []};
validator = {@(x) ischar(x) || iscellstr(x) || isnumeric(x), ...  % itemsToPrint
    @(x) ischar(x) & any(strcmp(x, [{'sum', 'all', 'obj', 'con', 'var'}, otherPrintTypes])), ...  % typeToPrint
    @(x) (isstruct(x) && (isfield(x, 'status') || isfield(x, 'stat')) || isvector(x)), ...  % solToPrint
    @(x) ischar(x) | iscellstr(x), ...  % colNames
    @(x) ischar(x) | iscellstr(x), ...  % rowNames
    @ischar, ...  % fileToPrint
    ... %%%%% add here other model format to support
    @(x) ischar(x) && any(strcmpi(x, {'cobra', 'gurobi', 'cplex'})), ...  % modelFormat  
    };

% get all potentially supplied COBRA parameter names
problemTypes = {'LP'};

[funParams, ~, ~] = parseCobraVarargin(varargin, optArgin, defaultValues, validator, problemTypes, '', true);

[itemsToPrint, typeToPrint, solToPrint, colNames, rowNames, fileToPrint, modelFormat] = deal(funParams{:});

% further handling
if isnan(typeToPrint)
    % if `typeToPrint` is not defined
    if isscalar(itemsToPrint) && isnumeric(itemsToPrint) && isnan(itemsToPrint)
        % if `itemsToPrint` is also not define, print summary
        typeToPrint = 'sum';
    elseif ischar(itemsToPrint) && strcmp(itemsToPrint, 'all')
        % print everything if `itemsToPrint` is 'all'
        typeToPrint = 'all';
    elseif isnumeric(itemsToPrint)
        % assume it refers to constraints if input is an index vector
        typeToPrint = 'con';
        % else if itemsToPrint is a cell array of strings, it will be handled later to detect `typeToPrint`
    end
elseif isnan(itemsToPrint)
    % Default print all items if a certain type to print is specified.
    itemsToPrint = 'all';
end
if ischar(itemsToPrint) && ~strcmp(itemsToPrint, 'all')
    % if `itemsToPrint` is a string of item name, put it in a cell
    itemsToPrint = {itemsToPrint};
end
if iscellstr(itemsToPrint)
    % `typeToPrint` input is overriden and determined by what are found in itemsToPrint
    typeToPrint = [];
end
%% unify format
if isempty(fileToPrint) || isequal(fileToPrint, 1)
    fid = 1;
    fileToPrint = '';
else
    fid = fopen(fileToPrint,'a');
end

checkFormat = false;
%%%%% add here other model format to support and check
if ~((strcmpi(modelFormat, 'cplex') && isa(model, 'Cplex')) || ...
        (strcmpi(modelFormat, 'cobra') && isfield(model, 'A')) || ...
        (strcmpi(modelFormat, 'gurobi') && isfield(model, 'A')))
    checkFormat = true;
else
    modelFormat = lower(modelFormat);
end
if checkFormat
    if isa(model, 'Cplex')
        modelFormat = 'cplex';
    elseif isfield(model, 'A')
        % gurobi or Cobra optimization model
        if any(isfield(model, {'b', 'c', 'osense', 'csense'}))
            modelFormat = 'cobra';
            if ~all(isfield(model, {'b', 'c', 'osense', 'csense', 'lb', 'ub'}))
                model = fillMissingFields(model, 'cobra');
            end
        else
            % assume it is gurobi since gurobi has 'A' as the only required field
            modelFormat = 'gurobi';
            if ~all(isfield(model, {'obj', 'sense', 'rhs', 'lb', 'ub', 'vtype', 'modelsense'}))
                model = fillMissingFields(model, 'gurobi');
            end
        end
    elseif isfield(model, 'S')
        % Cobra metabolic model
        modelFormat = 'cobraGEM';
        model = fillMissingFields(model, 'cobraGEM');
    else
        error('The input model structure must be a Cplex object, Gurobi or Cobra model structure for the constraints to be printed.')
    end
end

matchName = getMatchName();
if ~isempty(solToPrint) && isvector(solToPrint) && ~isstruct(solToPrint)
    % solution input is a vector but not solution structure (e.g., manually created solution)
    solToPrint = struct(matchName.(modelFormat).full, solToPrint);
end

if strcmp(modelFormat, 'cplex')
    % get the solution and model structures from the cplex construct
    if isempty(solToPrint) && isprop(model, 'Solution')
        solToPrint = model.Solution;
    end
    model = model.Model;
    if ~isempty(solToPrint)
        solToPrint.slack = [solToPrint.ax - model.lhs, model.rhs - solToPrint.ax];
    end
end

% check colNames and rowNames
%%%%% add here how to get the dimension and names for variables and constraints
switch modelFormat
    case 'cplex'
        if isempty(colNames) && isfield(model,'colname')
            colNames = cellstr(model.colname);
        end
        if isfield(model,'rowname')
            rowNames = cellstr(model.rowname);
        end
        [m, n] = size(model.A);
    case 'gurobi'
        if isempty(colNames) && isfield(model,'varnames')
            colNames = model.varnames(:);
        end
        if isfield(model, 'constrnames')
            rowNames = model.constrnames(:);
        end
        [m, n] = size(model.A);
    case 'cobraGEM'
        [m, n] = size(model.A);
        if isempty(colNames) && isfield(model, 'rxns')
            colNames = model.rxns;
        end
        if isempty(rowNames) && isfield(model, 'mets')
            rowNames = model.mets;
            if numel(rowNames) < m
                rowNames((end + 1):m) = strcat('Add_constr_', cellstr(num2str((1:(m - numel(rowNames)))')));
            end
        end
        % converting solution's fields
        if ~isempty(solToPrint)
            solToPrint.obj = solToPrint.f;
            solToPrint.full = solToPrint.v;
            if isfield(solToPrint, 'ctrs_y')
                solToPrint.dual = [solToPrint.y; solToPrint.ctrs_y];
            else
                solToPrint.dual = solToPrint.y;
            end
            if isfield(solToPrint, 'ctrs_slack')
                solToPrint.slack = [solToPrint.s; solToPrint.ctrs_slack];
            else
                solToPrint.slack = solToPrint.s;
            end
            solToPrint.rcost = solToPrint.w;
        end
        
        modelFormat = 'cobra';
    otherwise
        [m, n] = size(model.A);
end
if isempty(colNames)
    colNames = strcat('x',cellstr(num2str((1:n)','%-d')));
end
if isempty(rowNames)
    rowNames = strcat('r',cellstr(num2str((1:m)','%-d')));
end

%%%%% Replace here the symbols for the sign of constraints if not using
%%%%% 'E', 'L' and 'G' for equal to, less than or equal and greater than or
%%%%% equal to, respectively
if strcmp(modelFormat, 'gurobi') && isfield(model, 'sense')
    % replace the constraint sense symbols with lettters
    model.sense = strrep(strrep(strrep(model.sense(:)', '=', 'E'), '<', 'L'), '>', 'G')';
end

if ~isfield(model, matchName.(modelFormat).vartype.fieldname)
    % assume all variables are continuous if not specified
    model.(matchName.(modelFormat).vartype.fieldname) = char(matchName.(modelFormat).vartype.C * ones(n, 1));
end

%% to print certain properties or not
% properties of variables and constraints in solutions to be printed.
% Must be subset of properties in the `matchName` above
%%%%%%% Modify to add other properties
varProp = {'full', 'rcost'};
conProp = {'dual', 'slack'};

printVarProp = struct();
printConProp = struct();
for i = 1:numel(varProp)
    printVarProp.(varProp{i}) = false;
    if ~isempty(solToPrint)
        if isfield(solToPrint, matchName.(modelFormat).(varProp{i})) ...
                && numel(solToPrint.(matchName.(modelFormat).(varProp{i}))) == n
            printVarProp.(varProp{i}) = true;
        end
    end
end
for i = 1:numel(conProp)
    printConProp.(conProp{i}) = false;
    if ~isempty(solToPrint)
        if isfield(solToPrint, matchName.(modelFormat).(conProp{i})) ...
                && numel(solToPrint.(matchName.(modelFormat).(conProp{i}))) == m
            printConProp.(conProp{i}) = true;
        end
    end
end

colNames0 = colNames;
if printVarProp.full
    % append variable values in the solution to variable names
    colNames = strcat(colNames(:), ' (', numToFormattedString(solToPrint.(matchName.(modelFormat).full)), ')');
end

%% get sizes of other print types
%%%%%%%  add here the field names for other special constraints specific to certain model formats
nSC = struct();
otherPrintTypeKeys = struct();
otherPrintTypeKeys.sos = struct(...
    'gurobi', 'sos', ...
    'cplex', 'sos');
otherPrintTypeKeys.ind = struct(...
    'gurobi', 'genconind', ...
    'cplex', 'indicator');

for i = 1:numel(otherPrintTypes)
    if any(strcmp(typeToPrint, {'all', 'sum', otherPrintTypes{i}}))
        nSC.(otherPrintTypes{i}) = 0;
        if isfield(otherPrintTypeKeys.(otherPrintTypes{i}), modelFormat) && ...
                isfield(model, otherPrintTypeKeys.(otherPrintTypes{i}).(modelFormat))
            nSC.(otherPrintTypes{i}) = numel(model.(otherPrintTypeKeys.(otherPrintTypes{i}).(modelFormat)));
        end
    end
    if strcmp(typeToPrint, otherPrintTypes{i}) && nSC.(otherPrintTypes{i}) == 0
        warning('There are not any %s detected in the model.\n', otherPrintTypesNames.(otherPrintTypes{i}));
        return
    end
end

%% handle input itemsToPrint if it is a cell array of strings
if iscellstr(itemsToPrint)
    % find whether the input cells belong to the following categories
    allFound = false(numel(itemsToPrint), 1);
    [isCon, idCon] = ismember(itemsToPrint, rowNames);
    [isVar, idVar] = ismember(itemsToPrint, colNames0);
    isVar = isVar & ~isCon;
    idVar(~isVar) = 0;
    allFound = allFound | isCon | isVar;
    [isOtherItems, idOtherItems] = deal(struct());
    for i = 1:numel(otherPrintTypes)
        % search also other types if there are name fields in other items
        if isfield(otherPrintTypeKeys.(otherPrintTypes{i}), modelFormat) && ...
                isfield(model, otherPrintTypeKeys.(otherPrintTypes{i}).(modelFormat)) && ...
                isfield(model.(otherPrintTypeKeys.(otherPrintTypes{i}).(modelFormat)), 'name')
            [isOtherItems.(otherPrintTypes{i}), idOtherItems.(otherPrintTypes{i})] ...
                = ismember(itemsToPrint, {model.(otherPrintTypeKeys.(otherPrintTypes{i}).(modelFormat)).name});
            isOtherItems.(otherPrintTypes{i}) = isOtherItems.(otherPrintTypes{i}) & ~allFound;
            idOtherItems.(otherPrintTypes{i})(~isOtherItems.(otherPrintTypes{i})) = 0;
            allFound = allFound | isOtherItems.(otherPrintTypes{i});
        end
    end
    if ~all(allFound)
        error('The following items cannot be found in the model:\n%s\n', ...
            strjoin(itemsToPrint(~allFound), '\n'));
    end
    if ~isempty(fileToPrint)
        fclose(fid);
    end
    % print each item type in turn
    if any(isCon)
        printModel(model, idCon(isCon), 'con', solToPrint, colNames0, rowNames, fileToPrint, modelFormat);
    end
    for i = 1:numel(otherPrintTypes)
        if isfield(isOtherItems, otherPrintTypes{i}) && any(isOtherItems.(otherPrintTypes{i}))
            printModel(model, idOtherItems.(otherPrintTypes{i})(isOtherItems.(otherPrintTypes{i})), ...
                otherPrintTypes{i}, solToPrint, colNames0, rowNames, fileToPrint, modelFormat);
        end
    end
    if any(isVar)
        printModel(model, idVar(isVar), 'var', solToPrint, colNames0, rowNames, fileToPrint, modelFormat);
    end
    return
end

%% actual printing
switch typeToPrint
    case 'all'
        % print the entire model
        fprintf(fid, '%s model:\n', modelFormat);
        if ~isempty(fileToPrint)
          fclose(fid);
        end
        if strcmp(modelFormat, 'cplex')
            tmp = model;
            model = Cplex();
            model.Model = tmp;
            clear tmp
        end
        printModel(model, [], 'obj', solToPrint, colNames0, rowNames, fileToPrint, modelFormat)
        printModel(model, 1:m, 'con', solToPrint, colNames0, rowNames, fileToPrint, modelFormat)
        for i = 1:numel(otherPrintTypes)
            if nSC.(otherPrintTypes{i}) > 0
                printModel(model, 1:nSC.(otherPrintTypes{i}), otherPrintTypes{i}, ...
                    solToPrint, colNames0, rowNames, fileToPrint, modelFormat)
            end
        end
        printModel(model, 1:n, 'var', solToPrint, colNames0, rowNames, fileToPrint, modelFormat)
        return
        
    case 'sum'
        % print a summary only
        fprintf(fid, '%s model:\n', modelFormat);
        if ~isempty(fileToPrint)
            fclose(fid);
        end
        % print objective function
        printModel(model, [], 'obj', solToPrint, colNames0, rowNames, fileToPrint, modelFormat)
        if isempty(fileToPrint)
            fid = 1;
        else
            fid = fopen(fileToPrint,'a');
        end
        fprintf(fid, 'subject to \n   %d linear constraints\n', m);
        for i = 1:numel(otherPrintTypes)
            if nSC.(otherPrintTypes{i}) > 0
                fprintf(fid, '   %d %s constraints\n', nSC.(otherPrintTypes{i}), otherPrintTypesNames.(otherPrintTypes{i}));
            end
        end
        fprintf(fid, '   %d continuous, %d binary and %d general integer variables\n', ...
            sum(model.(matchName.(modelFormat).vartype.fieldname) == matchName.(modelFormat).vartype.C), ...
            sum(model.(matchName.(modelFormat).vartype.fieldname) == matchName.(modelFormat).vartype.B), ...
            sum(model.(matchName.(modelFormat).vartype.fieldname) == matchName.(modelFormat).vartype.I));
        if ~isempty(fileToPrint)
            fclose(fid);
        end
        return
        
    case 'obj'
        % print objective function
        if isequal(model.(matchName.(modelFormat).maxmin.fieldname), matchName.(modelFormat).maxmin.max)
            s = 'max';
        elseif isequal(model.(matchName.(modelFormat).maxmin.fieldname), matchName.(modelFormat).maxmin.min)
            s = 'min';
        else
            error('The optimization sense (model.%s) should be either %s or %s for a %s model.', ...
                matchName.(modelFormat).maxmin.fieldname, matchName.(modelFormat).maxmin.min, matchName.(modelFormat).maxmin.max, modelFormat)
        end
        objVector = model.(matchName.(modelFormat).obj);
        fprintf(fid, '%s  %s', s, printLinearCombination(colNames, find(objVector), objVector(objVector ~= 0)));
        if ~isempty(solToPrint)
            if isfield(solToPrint, matchName.(modelFormat).f) && ~isnan(solToPrint.(matchName.(modelFormat).f))
                fprintf(fid, '  (= %s)', numToFormattedString(solToPrint.(matchName.(modelFormat).f)));
            else
                fprintf(fid, '  (= NaN)');
            end
        end
        fprintf('\n');
                
    case 'con'
        % print constraints
        if ischar(itemsToPrint) && strcmp(itemsToPrint, 'all')
            % print all constraints
            itemsToPrint = 1:m;
        elseif any(itemsToPrint ~= floor(itemsToPrint) | itemsToPrint < 1 | itemsToPrint > m)
            % check index vector input  
            error('Constraint index vector must contain only integers between 1 and the number of constraints (%d detected)', m)
        end
        printBothSides = isfield(matchName.(modelFormat), 'lhs');
        csenseArrows = struct('E', '=', ...
            'L', '<=', ...
            'G', '>=');
        for i = itemsToPrint(:)'
            fprintf(fid, 'con %d %s:  ', i, rowNames{i});
            % print LHS
            if printBothSides
                fprintf(fid, '%s <= ', numToFormattedString(model.(matchName.(modelFormat).lhs)(i)));
            end
            % print the linear constraint
            a = model.(matchName.(modelFormat).A)(i, :);
            fprintf(fid, '%s ', printLinearCombination(colNames, find(a), a(a ~= 0)));
            % print RHS
            if printBothSides
                s = '<=';
            else
                s = csenseArrows.(model.(matchName.(modelFormat).csense)(i));
            end
            fprintf(fid, '%s %s', s, numToFormattedString(model.(matchName.(modelFormat).rhs)(i)));
            if printConProp.dual || printConProp.slack
                fprintf(fid, ' (');
                if printConProp.dual
                    fprintf(fid, 'dual=%s', numToFormattedString(solToPrint.(matchName.(modelFormat).dual)(i)));
                    if printConProp.slack
                        fprintf(fid, ', ');
                    end
                end
                if printConProp.slack
                    fprintf(fid, 'slack=%s', strjoin(numToFormattedString(solToPrint.(matchName.(modelFormat).slack)(i, :), [], [], [], true), ','));
                end
                fprintf(fid, ')');
            end
            fprintf(fid, '\n');
        end

    case 'var'
        % print variables
        if ischar(itemsToPrint) && strcmp(itemsToPrint, 'all')
            itemsToPrint = 1:n;
        elseif any(itemsToPrint ~= floor(itemsToPrint) | itemsToPrint < 1 | itemsToPrint > n)
            error('Variable index vector must contain only integers between 1 and the number of variables (%d detected)', n)
        end
        varTypeNames = struct(matchName.(modelFormat).vartype.C, 'continuous', ...
            matchName.(modelFormat).vartype.B, 'binary', ...
            matchName.(modelFormat).vartype.I, 'integer');
        for i = itemsToPrint(:)'
            fprintf(fid, 'var %d, %s: %s <= %s <= %s', i, varTypeNames.(model.(matchName.(modelFormat).vartype.fieldname)(i)), ...
                numToFormattedString(model.lb(i)), colNames{i}, numToFormattedString(model.ub(i)));
            if printVarProp.rcost
                fprintf(fid, ' (rcost=%s)', numToFormattedString(solToPrint.(matchName.(modelFormat).rcost)(i)));
            end
            fprintf(fid, '\n');
        end
    otherwise
        % print other items
        if strcmp(itemsToPrint, 'all')
            itemsToPrint = 1:nSC.(typeToPrint);
        end
        if any(itemsToPrint ~= floor(itemsToPrint) | itemsToPrint < 1 | itemsToPrint > nSC.(typeToPrint))
            error('%s index vector must contain only integers between 1 and the number of entries (%d detected)', nSC.(typeToPrint))
        end
        %%%%%%%%%%%  add here the printing methods to new special items
        %%%%%%%%%%%  `printFields` are the fields in the structure for the
        %%%%%%%%%%%  special constraint that will be printed, %i is a
        %%%%%%%%%%%  reserved word for the index of that constraints
        %%%%%%%%%%%  `printFun` is a function handle to print the constraints,
        %%%%%%%%%%%  where the input is a cell array of fields in printFields
        %%%%%%%%%%%  in the structure for the special constraintswith. See the examples below. 
        switch [modelFormat, ',', typeToPrint]
            case 'gurobi,sos'
                printFields = {'%i', 'index', 'type'};
                printFun = @(x) fprintf(fid, 'sos %d (type %d):  %s\n', ...
                    x{1}, x{3}, strjoin(colNames(x{2}), ', '));
            case 'gurobi,ind'
                printFields = {'%i', 'binvar', 'binval', 'a', 'sense', 'rhs'};
                printFun = @(x) fprintf(fid, 'Ind %d:  %s = %d  =>  %s %s %s\n', ...
                    x{1}, colNames{x{2}}, x{3}, ...
                    printLinearCombination(colNames, find(x{4}), x{4}(x{4} ~= 0)), ...
                    x{5}, numToFormattedString(x{6}));
            case 'cplex,sos'
                printFields = {'%i', 'name', 'var', 'type'};
                printFun = ...
                    @(x) fprintf(fid, 'sos %d %s (type %d):  %s\n', x{1}, x{2}, x{4}, strjoin(colNames(x{3}), ', '));
            case 'cplex,ind'
                cplexIndType = {'=>', '<=', '<=>'};
                printFields = {'name', 'variable', 'complemented', 'type', 'a', 'sense', 'rhs', '%i'};
                printFun = @(x) fprintf(fid, 'ind %d %s:  %s = %d  %s  %s %s %s\n', ...
                    x{8}, x{1}, colNames{x{2}}, 1 - x{3}, cplexIndType{x{4}}, ...
                    printLinearCombination(colNames, find(x{5}), x{5}(x{5} ~= 0)), ...
                    x{6}, numToFormattedString(x{7}));
            otherwise
                error('Printing %s is currently not supported for %s model input.', otherPrintTypesNames.(typeToPrint))
        end
        for i = itemsToPrint(:)'
            s = model.(otherPrintTypeKeys.(typeToPrint).(modelFormat))(i);
            x = cell(numel(printFields), 1);
            for j = 1:numel(printFields)
                if strcmp(printFields{j}, '%i')
                    x{j} = i;
                else
                    x{j} = s.(printFields{j});
                end
            end
            printFun(x);
        end

end
if ~isempty(fileToPrint)
    fclose(fid);
end

end

function matchName = getMatchName()
%% match field names used for printing linear constraints and variables
%%%%%%  add here to match the used field names for other model format
matchName = struct();
% constraint matrix
matchName.gurobi.A = 'A';
matchName.cplex.A = 'A';
matchName.cobra.A = 'A';
% RHS
matchName.gurobi.rhs = 'rhs';
matchName.cplex.rhs = 'rhs';
matchName.cobra.rhs = 'b';
% LHS (Cplex only for now)
matchName.cplex.lhs = 'lhs';
% objective function
matchName.gurobi.obj = 'obj';
matchName.cplex.obj = 'obj';
matchName.cobra.obj = 'c';
% constraint sense (not needed if the model is like Cplex using LHS <= Ax <= RHS)
matchName.gurobi.csense = 'sense';
matchName.cobra.csense = 'csense';

% variable type
matchName.gurobi.vartype.fieldname = 'vtype';
matchName.gurobi.vartype.C = 'C';
matchName.gurobi.vartype.B = 'B';
matchName.gurobi.vartype.I = 'I';

matchName.cobra.vartype.fieldname = 'vartype';
matchName.cobra.vartype.C = 'C';
matchName.cobra.vartype.B = 'B';
matchName.cobra.vartype.I = 'I';

matchName.cplex.vartype.fieldname = 'ctype';
matchName.cplex.vartype.C = 'C';
matchName.cplex.vartype.B = 'B';
matchName.cplex.vartype.I = 'I';

% maximization or minimization
matchName.gurobi.maxmin.fieldname = 'modelsense';
matchName.gurobi.maxmin.max = 'max';
matchName.gurobi.maxmin.min = 'min';

matchName.cplex.maxmin.fieldname = 'sense';
matchName.cplex.maxmin.max = 'maximize';
matchName.cplex.maxmin.min = 'minimize';

matchName.cobra.maxmin.fieldname = 'osense';
matchName.cobra.maxmin.max = -1;
matchName.cobra.maxmin.min = 1;

% solution data
% objective function value
matchName.gurobi.f = 'objval';
matchName.cplex.f = 'objval';
matchName.cobra.f = 'obj';
% full solution vector
matchName.gurobi.full = 'x';
matchName.cplex.full = 'x';
matchName.cobra.full = 'full';
% reduced cost
matchName.gurobi.rcost = 'rc';
matchName.cplex.rcost = 'reducedcost';
matchName.cobra.rcost = 'rcost';
% dual values / shadow price
matchName.gurobi.dual = 'pi';
matchName.cplex.dual = 'dual';
matchName.cobra.dual = 'dual';
% slackness of constraints (Ax - b for >= constraints and b - Ax for <= constraints)
matchName.gurobi.slack = 'slack';
matchName.cplex.slack = 'slack';
matchName.cobra.slack = 'slack';
end

function model = fillMissingFields(model, modelFormat)
% Add missing fields in the model required for printing the model
% add fields for consistency with `modelFormat` = `cobra`
if strcmp(modelFormat, 'cobraGEM')
    if isfield(model, 'C')
        % additional constraints in a cobra metabolic model
        model.A = [model.S; model.C];
        if ~isfield(model, 'csense')
            model.csense = char('E' * ones(size(model.S, 1), 1));
        end
        if ~isfield(model, 'b')
            model.b = zeros(size(model.S, 1), 1);
        end
        if isfield(model, 'd')
            model.b = [model.b; model.d];
        else
            model.b = [model.b; zeros(size(model.C, 1), 1)];
        end
        if isfield(model, 'dsense')
            model.csense = [model.csense(:); model.dsense(:)];
        else
            model.csense = [model.csense(:); char('E' * ones(size(model.C, 1), 1))];
        end
    else
        model.A = model.S;
    end
    if isfield(model, 'osenseStr')
        if strcmp(model.osenseStr, 'max')
            model.osense = -1;
        elseif strcmp(model.osenseStr, 'min')
            model.osense = 1;
        else
            error('Cobra metabolic model should have .osenseStr equal to either ''max'' or ''min''.')
        end
    end
end
%%%%%%%% Add here the field names and default values for objective
%%%%%%%% function, constraint sense, RHS, LB, UB, variable type and
%%%%%%%% optimization sense for other model formats.
switch modelFormat
    case 'gurobi'
        [m, n] = size(model.A);
        % default values based on Gurobi manual
        defaultValues = {'obj', zeros(n, 1); ...
            'sense', char('<' * ones(m, 1)); ...
            'rhs', zeros(m, 1); ...
            'lb', zeros(n, 1); ...
            'ub', inf(n, 1); ...
            'vtype', char('C' * ones(n, 1)); ...
            'modelsense', 'min'};
    case {'cobra', 'cobraGEM'}
        [m, n] = size(model.A);
        % default values based on solveCobraLP
        defaultValues = {'c', NaN(n, 1); ...
            'csense', char('E' * ones(m, 1)); ...
            'b', zeros(m, 1); ...
            'lb', NaN(n, 1); ...
            'ub', NaN(n, 1); ...
            'vartype', char('C' * ones(n, 1)); ...
            'osense', -1};
end
for i = 1:size(defaultValues, 1)
    if ~isfield(model, defaultValues{i, 1})
        model.(defaultValues{i, 1}) = defaultValues{i, 2};
    end
end
end

function s = printLinearCombination(names, index, coeff)
% generate a string concatenating a cell array of strings or linear
% combination of them.
%
% USAGE:
%    s = printLinearCombination(names, index, coeff)
%
% INPUTS:
%    names:      cell array of strings to be concatenated and printed
%    index:      index vector for names to be extracted. `names(index)` will be concatenated
%    coeff:      coefficent vector with the same length as `index` for
%                generating the linear combination text
%
% OUTPUT:
%    s:          the corresponding generated string

if ~isempty(index) && all(isnan(coeff))
    s = 'NaN';
    return
end
coeff = full(coeff);    
s = '';
for jC = 1:numel(index)
    sJ = numToFormattedString(coeff(jC));
    % replace the resultant string by '' if the result is '1', or by '-' if it is '-1'.
    if strcmp(sJ, '1')
        sJ = '';
    elseif strcmp(sJ, '-1')
        sJ = '-';
    end
    % if it is not the first term in the linear combination, add space and/or sign
    if jC > 1
        if coeff(jC) < 0
            sJ = strrep(sJ, '-', ' - ');
        else
            sJ = [' + ' sJ];
        end
    end

    s = [s, sJ, ' ', strtrim(names{index(jC)})];
end
if isempty(s)
    % if everything is zero
    s = '0';
end
end
