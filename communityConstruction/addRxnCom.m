function [modelCom, id] = addRxnCom(modelCom,rxnID, rxnName, metaboliteList,...
    stoichCoeffList,revFlag,subSystem,grRule,checkDuplicate, addMets)
%[modelCom, id] = addRxnCom(modelCom,rxnID, rxnName, metaboliteList,...
%     stoichCoeffList,revFlag,subSystem,grRule,checkDuplicate, addMets)
% Add a reaction into the community model. It intends to have the similar
% input format as addReaction
%
%Input:
%   modelCom:       community model
%   rxnID:          name in modelCom.rxns
%   rxnName:        name in modelCom.rxnNames
%   metaboliteList: a reaction formula or list of metabolites in the
%                   reaction. If it is a reaction formula, later arguments
%                   are optional or stoichCoeffList must be left empty
%   stoichCoeffList:stoichiometry coefficient array corresponding to metaboliteList 
%   revFlag:        true for reversible, false for irreversible
%   subSystem:      name in modelCom.subSystem
%   grRule:         rule in modelCom.grRule
%   checkDuplicate: check if the same reaction exists. Default true
%   addMets:        automatically add new metabolites in the reaction using
%                   addMetCom. Default true. Error if false and new mets
%                   are found.
%
%Siu Hung Joshua Chan Nov 2016

parseForm = false;
if nargin < 5
    %input formula
    parseForm = true;
elseif isempty(stoichCoeffList)
    parseForm = true;
end
if parseForm
    [metaboliteList,stoichCoeffList,revFlag0] = parseRxnFormula(metaboliteList);
    if ~exist('revFlag','var')
        revFlag = revFlag0;
    elseif isempty(revFlag)
        revFlag = revFlag0;
    end
end
if numel(metaboliteList) ~= numel(stoichCoeffList)
    error('Sizes of metaboliteList and stoichCoeffList not equal.');
end
if ~exist('subSystem','var')
    subSystem = '';
elseif isempty(subSystem)
    subSystem = '';
end
if ~exist('grRule','var')
    grRule = '';
elseif isempty(grRule)
    grRule = '';
end
if ~exist('checkDuplicate','var')
    checkDuplicate = true;
elseif isempty(checkDuplicate)
    checkDuplicate = true;
end
if ~exist('addMets','var')
    addMets = true;
elseif isempty(addMets)
    addMets = true;
end

if isnumeric(metaboliteList)
    metIDs = metaboliteList;
else
    metIDs = findMetIDs(modelCom,metaboliteList);
end
if ~all(metIDs)
    if addMets
        met2add = find(metIDs == 0);
        for j = 1:numel(met2add)
            modelCom = addMetCom(modelCom, metaboliteList{met2add(j)});
        end
    else
        error('Metabolites in the reaction formula not in the model.')
    end
    metIDs(met2add) = findMetIDs(modelCom, metaboliteList(met2add));
end
speciesId = unique(modelCom.indCom.metSps(metIDs));

[m,n] =size(modelCom.S);
nSp = numel(modelCom.infoCom.spAbbr);
%community exchange reaction or not
comEx = false;
%species exchange reaction or not
spEx = false;
%assign Id
if numel(speciesId) > 1
    if numel(speciesId) == 2 && (speciesId(1) == 0 || speciesId(2) == 0)
    %exchange between a species and the community space, species specific
        speciesId = speciesId(speciesId ~= 0);
        if numel(stoichCoeffList) == 2
            %exchange reaction
            spEx = true;
        end
    else
    %interspecies reaction. Give speciesId -1. Add before community reactions
        speciesId = -1;
        id = find(modelCom.indCom.rxnSps == 0, 1);
    end
end

if speciesId >= 0 %species-specific or community reaction
    if speciesId == 0 %community reaction
        if numel(metIDs) == 1  %community exchange
            comEx = true;
            if stoichCoeffList < 0 % export
                id = n + 1; %put at the end
            elseif stoichCoeffList > 0 % uptake
                id = modelCom.indCom.EXcom(1,2); %put just before export
            end
            if isempty(subSystem)
                subSystem = 'Community exchange';
            end
        else
            id = modelCom.indCom.EXcom(1,1); %put just before uptake
        end
    elseif speciesId == nSp       
        %reaction within the last species
        %put before first interspecies or community reaction
        id = find(modelCom.indCom.rxnSps <= 0, 1); 
    else
        %reaction within other single species
        id = find(modelCom.indCom.rxnSps == speciesId + 1, 1);
    end
end

%check rxnID
%check any species's name attached to the end of rxnID
spId = find(cellfun(@(x) regexpcmp(rxnID, [x '$']), modelCom.infoCom.spAbbr));
    %spId = find(strcmp(modelCom.infoCom.spAbbr, spId),1);
addSpeciesTag = isempty(spId);
        
%check if it is designated as a community reaction
if addSpeciesTag && length(rxnID) > 3
    if strcmp(rxnID(end-2:end), '(u)')
        spId = 0;
        addSpeciesTag = false;
    end
end
%add species identifier
if addSpeciesTag
    if speciesId == -1
        speciesTag = '';
    elseif speciesId == 0
        speciesTag = '(u)';
    else
        speciesTag = ['_' modelCom.infoCom.spAbbr{speciesId}];
    end
    rxnID = [rxnID speciesTag];
else
    if speciesId ~= spId
        error('SpeciesId inferred from metabolites = %d. From rxnID = %d. Not the same.',...
            speciesId, spId);
    end
end

if any(strcmp(modelCom.rxns, rxnID))
    error('rxnID %s already exists in the model.', rxnID);
end
newCol = sparse(metIDs, ones(numel(metIDs),1), stoichCoeffList, m, 1);
if checkDuplicate
    yn = ismember(modelCom.S(newCol ~= 0,:)', newCol(newCol ~= 0,:)','rows');
    if any(yn)
        f = ~any(modelCom.S(newCol == 0, yn),1);
        if any(f)
            ind = find(yn);
            f = ind(f);
            errorMsg = sprintf('The reaction is the same as ');
            for j = 1:(numel(f)-1)
                errorMsg = [errorMsg sprintf('%s, ', modelCom.rxns{f(j)})];
            end
            errorMsg = [errorMsg sprintf('%s.\n', modelCom.rxns{f(end)})];
            error(errorMsg);
        end
    end
end

%update S matrix
modelCom.S(:,(id + 1): (n + 1)) = modelCom.S(:,id:n);
modelCom.S(:,id) = 0;
modelCom.S(metIDs,id) = stoichCoeffList(:);
%update main fields
field2change = {'rxns'; 'rxnNames'; 'rev'; 'c'; 'lb'; 'ub'; 'subSystems'; 'grRules'};
fieldValue = {rxnID; rxnName; revFlag; 0; revFlag * -1000; 1000; subSystem; grRule};
for j = 1:numel(field2change)
    if isfield(modelCom,field2change{j})
        modelCom.(field2change{j})((id + 1): (n + 1)) = modelCom.(field2change{j})(id:n);
        if iscell(modelCom.(field2change{j}))
            modelCom.(field2change{j}){id} = fieldValue{j};
        else
            modelCom.(field2change{j})(id) = fieldValue{j};
        end
    end
end
modelCom.indCom.rxnSps((id + 1): (n + 1)) = modelCom.indCom.rxnSps(id:n);
modelCom.indCom.rxnSps(id) = speciesId;
%update reaction indices, add 1 to all reactions with indices >= id 
modelCom.indCom.EXcom(modelCom.indCom.EXcom >= id) = modelCom.indCom.EXcom(modelCom.indCom.EXcom >= id) + 1;
modelCom.indCom.EXsp(modelCom.indCom.EXsp >= id) = modelCom.indCom.EXsp(modelCom.indCom.EXsp >= id) + 1;
modelCom.indCom.spBm(modelCom.indCom.spBm >= id) = modelCom.indCom.spBm(modelCom.indCom.spBm >= id) + 1;
if comEx
    metCom = modelCom.indCom.metSps == 0;
    metComId = sum(metCom(1:metIDs));
    if stoichCoeffList < 0
        if modelCom.indCom.EXcom(metComId, 2) > 0
            metCom = find(metCom);
            error('The community export reaction for %s already exists as %s.',...
                modelCom.mets{metCom(metComId)}, modelCom.rxns{modelCom.indCom.EXcom(metComId, 2)});
        else
            modelCom.indCom.EXcom(metComId, 2) = id;
        end
    else
        if modelCom.indCom.EXcom(metComId, 1) > 0
            metCom = find(metCom);
            error('The community uptake reaction for %s already exists as %s.',...
                modelCom.mets{metCom(metComId)}, modelCom.rxns{modelCom.indCom.EXcom(metComId, 2)});
        else
            modelCom.indCom.EXcom(metComId, 1) = id;
        end
    end
end
if spEx
    speciesId = modelCom.indCom.metSps(metIDs);
    metCom = modelCom.indCom.metSps == 0;
    metComId = sum(metCom(1:metIDs(speciesId == 0)));
    if modelCom.indCom.EXsp(metComId, speciesId(speciesId ~= 0)) > 0
        metCom = find(metCom);
        error('The species-community exchange reaction for %s already exists as %s.',...
            strrep(modelCom.mets{metCom(metComId)},'[u]',''),...
            modelCom.rxns{modelCom.indCom.EXsp(metComId, speciesId(speciesId ~= 0))});
    else
        modelCom.indCom.EXsp(metComId, speciesId(speciesId ~= 0)) = id;
    end
end


%move the position of rxn-related fields
%the new value needs to be filled in manually
field = setdiff(fieldnames(modelCom),field2change);
for j = 1:numel(field)
    if length(field{j}) >= 3 
        if strcmp(field{j}(1:3), 'rxn')
            if size(modelCom.(field{j}), 1) == n
                modelCom.(field{j})((id + 1): (n + 1),:) = modelCom.(field{j})(id:n,:);
                if iscell(modelCom.(field{j}))
                    modelCom.(field{j})(id,:) = repmat({''}, 1, size(modelCom.(field{j}),2));
                else
                    modelCom.(field{j})(id,:) = 0;
                end
            end
            if size(modelCom.(field{j}), 2) == n
                modelCom.(field{j})(:,(id + 1): (n + 1)) = modelCom.(field{j})(:,id:n);
                if iscell(modelCom.(field{j}))
                    modelCom.(field{j})(:,id) = repmat({''}, size(modelCom.(field{j}),1), 1);
                else
                    modelCom.(field{j})(:,id) = 0;
                end
            end
        end
    end
end

% update .infoCom using .indCom
modelCom.infoCom = infoCom2indCom(modelCom, modelCom.indCom, true, modelCom.infoCom.spAbbr, modelCom.infoCom.spName);
end
    