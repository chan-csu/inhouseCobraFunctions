function varargout = isBlocked(model, metrxn, useModelExCond)
% Test whether a metabolite or reaction is blocked.
%
% USAGE:
%    [isPro, isCon] = isBlocked(model, mets, useModelExCond)
%    [isFor, isRev] = isBlocked(model, rxns, useModelExCond)
%
% INPUT:
%    model:          COBRA model
%
% OPTIONAL INPUTS:
%    mets/rxns:      cell array of metabolite/reaction IDs (default model.mets)
%    useModelExCond: 1 to use the uptake/export condition in the model (defined by model.lb and model.ub, default)
%                    0 to allow uptake and export of anything from the model, the most relaxed situation
% 
% OUTPUTS:
%    isPro:          true if the metabolite is producible by the model (an added sink reaction can have positive flux)  
%    isCon:          true if the metabolite is consumable by the model (an added sink reaction can have negative flux)  
%    isFor:          true if the reaction can carry positive flux
%    isCon:          true if the reaction can carry negative flux

if nargin < 2 || isempty(metrxn)
    metrxn = model.mets;
end
if nargin < 3 || isempty(useModelExCond)
    useModelExCond = 1;
end
if ischar(metrxn)
    metrxn = {metrxn};
end

% check whether metrxn is mets or rxns
isMet = all(findMetIDs(model, metrxn));
if ~isMet && ~all(findRxnIDs(model, metrxn))
    error('the 2nd argument `metrxn` must be a cell array of either all metabolite IDs or all reaction IDs')
end

if ~useModelExCond
    % relax the uptake and export conditions
    rxnEx = sum(model.S ~= 0, 1) == 1;
    model.lb(rxnEx) = -1000;
    model.ub(rxnEx) = 1000;
end

% objective function not needed
model = changeObjective(model, model.rxns(1), 0);

varargout = {false(numel(metrxn), 1), false(numel(metrxn), 1)};
if isMet
    for j = 1:numel(metrxn)
        % add sink reaction
        modelJ = addSinkReactions(model, metrxn(j));
        % constrain positive flux
        modelJ = changeRxnBounds(modelJ, ['sink_' metrxn{j}],  0.1, 'l');
        s = optimizeCbModel(modelJ);
        if s.stat == 1  % if feasible
            varargout{1}(j) = true;
        end
        % constrain negative flux
        modelJ = changeRxnBounds(modelJ, ['sink_' metrxn{j}],  -1000, 'l');
        modelJ = changeRxnBounds(modelJ, ['sink_' metrxn{j}],  -0.1, 'u');
        s = optimizeCbModel(modelJ);
        if s.stat == 1  % if feasible
            varargout{2}(j) = true;
        end
    end
else
    for j = 1:numel(metrxn)
        % constrain positive flux
         modelJ = changeRxnBounds(model, metrxn{j},  0.1, 'l');
         modelJ = changeRxnBounds(modelJ, metrxn{j},  1000, 'u');
         s = optimizeCbModel(modelJ);
         if s.stat == 1  % if feasible
             varargout{1}(j) = true;
         end
         % constrain negative flux
         modelJ = changeRxnBounds(model, metrxn{j},  -1000, 'l');
         modelJ = changeRxnBounds(modelJ, metrxn{j},  -0.1, 'u');
         s = optimizeCbModel(modelJ);
         if s.stat == 1  % if feasible
             varargout{2}(j) = true;
         end
    end
end

end