function [modelCom] = removeRxnCom(modelCom, rxnID, rxnFieldKnown)
%[modelCom] = removeRxnCom(modelCom, rxnID, rxnFieldKnown)
% Add a metabolite to the model.
% id is the index of the metabolite inserted
% Update the field mets, metNames, metFormulas, metSps, b, and other fields
% starting with 'met' and have at least one dimension the same as modelCom.mets
% Other fields need to be mannually updated using 'id'

% if isnumeric(rxnID) && numel(rxnID) > 1
%     if isnumeric(rxnID)
%         rxnID = modelCom.rxns(rxnID);
%     end
%     for j = 1:numel(rxnID)
%         modelCom = removeRxnCom(modelCom,rxnID{j});
%     end
%     return
% end
if nargin < 3
    rxnFieldKnown = {};
end
if ~isnumeric(rxnID)
    rxnID = findRxnIDs(modelCom,rxnID);
    if ~all(rxnID)
        warning('Some reaction IDs are not present in the model.');
        return
    end
end

[~,nRxns] = size(modelCom.S);
if isfield(modelCom,'indCom')
    indCom = modelCom.indCom;
else
    indCom = infoCom2indCom(modelCom);
end
%store the correct reaction names
spBm = modelCom.rxns(indCom.spBm);
EXcom = cell(size(indCom.EXcom));
EXcom(indCom.EXcom~=0) = modelCom.rxns(indCom.EXcom(indCom.EXcom~=0));
EXsp = cell(size(indCom.EXsp));
EXsp(indCom.EXsp~=0) = modelCom.rxns(indCom.EXsp(indCom.EXsp~=0));
if isfield(indCom,'spATPM')
    spATPM = modelCom.rxns(indCom.spATPM);
end
%remove the reaction if it is in the above three categories
% spBm(indCom.spBm == rxnID) = {[]};
% EXcom(indCom.EXcom == rxnID) = {[]};
% EXsp(indCom.EXsp == rxnID) = {[]};
%sink and demand


% Construct vector to select rxns to be included in the model rapidly
selectRxns = true(nRxns,1);
selectRxns(rxnID) = false;

modelCom.S = modelCom.S(:,selectRxns);
modelCom.rxns = modelCom.rxns(selectRxns);
modelCom.lb = modelCom.lb(selectRxns);
modelCom.ub = modelCom.ub(selectRxns);
modelCom.rev = modelCom.rev(selectRxns);
if (isfield(modelCom,'c'))
    modelCom.c = modelCom.c(selectRxns);
end
if isfield(modelCom,'rxnGeneMat')
    modelCom.rxnGeneMat = modelCom.rxnGeneMat(selectRxns,:);
end
field = fieldnames(modelCom);
rxnField = union(field(strncmp('rxn',field,3)),rxnFieldKnown);
reactionFields = {'subSystems', 'rxnNames', 'rxnReferences', 'rxnECNumbers',...
    'ecNumbers', 'rxnNotes', 'confidenceScores', 'citations', 'rxnKeggID', ...
    'comments','grRules','rules'};
rxnField = intersect(union(reactionFields,rxnField),field);
rxnField = setdiff(rxnField,{'rxns','lb','ub','rev','c','rxnGeneMat'});
for i=1:length(rxnField)
	modelCom.(rxnField{i}) = modelCom.(rxnField{i})(selectRxns);
end

%update reaction indices
% indLogic = ismember(indCom.spBm, rxnID);
% indCom.spBm(indLogic) = 0;
% indCom.spBm(~indLogic) = findRxnIDs(modelCom, spBm(~indLogic));
% indLogic = ismember(indCom.EXcom, rxnID) | indCom.EXcom == 0;
% indCom.EXcom(indLogic) = 0;
% indCom.EXcom(~indLogic) = findRxnIDs(modelCom,EXcom(~indLogic));
% indLogic = ismember(indCom.EXsp, rxnID) | indCom.EXsp == 0;
% indCom.EXsp(indLogic) = 0;
% indCom.EXsp(~indLogic) = findRxnIDs(modelCom,EXsp(~indLogic));
selectRxnIDs = find(selectRxns);
[~, indCom.spBm] = ismember(indCom.spBm, selectRxnIDs);
[~, indCom.EXcom] = ismember(indCom.EXcom, selectRxnIDs);
[~, indCom.EXsp] = ismember(indCom.EXsp, selectRxnIDs);
if isfield(indCom,'spATPM')
    [~, indCom.spATPM] = ismember(indCom.spATPM, selectRxnIDs);
    %     indLogic = ismember(indCom.spATPM, rxnID);
    %     indCom.spATPM(indLogic) = 0;
    %     indCom.spATPM(~indLogic) = findRxnIDs(modelCom, spATPM(~indLogic));
end
indCom.rxnSps(rxnID) = [];
if isfield(indCom,'rxnSD')
    indCom.rxnSD = sum(modelCom.S ~= 0, 1) == 1;
    indCom.rxnSD(indCom.EXcom(indCom.EXcom~=0)) = false;
    indCom.rxnSD = find(indCom.rxnSD);
end

modelCom.indCom = indCom;
modelCom.infoCom = infoCom2indCom(modelCom,indCom,true,modelCom.infoCom.spAbbr,modelCom.infoCom.spName);
end
                
            


