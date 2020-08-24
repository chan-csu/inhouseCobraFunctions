function [modelCom, id] = addMetCom(modelCom, metID, metName, formula, speciesId)
%[modelCom, id] = addMetCom(modelCom, metID, metName, formula, speciesId)
% Add a metabolite to the model.
% id is the index of the metabolite inserted
% Update the field mets, metNames, metFormulas, metSps, b, and other fields
% starting with 'met' and have at least one dimension the same as modelCom.mets
% Other fields need to be mannually updated using 'id'
%
%Siu Hung Joshua Chan Nov 2016

underscore = regexp(metID, '_');
addSpeciesIden = true;
spId = [];
%check any species's name attached to the end of metID
if ~isempty(underscore)
    if length(metID) > underscore(end)
        spId = metID((underscore(end) + 1):length(metID));
        spId = find(strcmp(modelCom.infoCom.spAbbr, spId),1);
    end
end
%check if it is a community metabolite
if isempty(spId)
    if length(metID) > 3
        if strcmp(metID(end-2:end), '[u]')
            spId = 0;
        end
    end
end
%if not found in metID, need to supply speciesID
if ~isempty(spId)
    addSpeciesIden = false;
end
if nargin < 5
    if addSpeciesIden
        error('Cannot identify species ID for the input metabolite.');
    else
        %use speciesID detected from metID
        speciesId = spId;
    end
else
    %get a numeric species Id
    if iscell(speciesId)
        speciesId = speciesId{1};    
    end
    if ischar(speciesId)
        speciesId = find(strcmp(modelCom.infoCom.spAbbr, speciesId),1);
    end
    if ~addSpeciesIden
        if spId ~= speciesId
            error('Invalid metID, coinciding with species''s name')
        end
    end
end
if nargin < 4
    formula = '';
elseif isempty(formula)
    formula = '';
end
if nargin < 3
    metName = '';
elseif isempty(metName)
    metName = '';
end
if addSpeciesIden
    metID = [metID '_' modelCom.infoCom.spAbbr{speciesId}];
end
if any(strcmp(modelCom.mets, metID))
    error('metID %s already exists in the model.', metID);
end

nSp = numel(modelCom.infoCom.spAbbr);
nCom = size(modelCom.indCom.EXcom,1);
if speciesId == 0
    %community metabolites put at the end
    id = numel(modelCom.mets) + 1;
    %add row to modelCom.EXcom and modelCom.EXsp for adding indicies
    modelCom.indCom.EXcom(nCom + 1,:) = 0;
    modelCom.indCom.EXsp(nCom + 1,:) = 0;
else
    %other species's metabolites appended to its own set of metabolites
    id = find(modelCom.indCom.metSps == mod(speciesId + 1, nSp + 1), 1);
end

[m, n] =size(modelCom.S);
modelCom.S((id + 1): (m + 1), :) = modelCom.S(id:m,:);
modelCom.S(id,:) = 0;
modelCom.mets((id + 1): (m + 1)) = modelCom.mets(id:m);
modelCom.mets{id} = metID;
if isfield(modelCom, 'metNames')
    modelCom.metNames((id + 1): (m + 1)) = modelCom.metNames(id:m);
    modelCom.metNames{id} = metName;
end
if isfield(modelCom, 'metFormulas')
    modelCom.metFormulas((id + 1): (m + 1)) = modelCom.metFormulas(id:m);
    modelCom.metFormulas{id} = formula;
end

if isfield(modelCom, 'b')
    modelCom.b((id + 1): (m + 1)) = modelCom.b(id:m);
    modelCom.b(id) = 0;
end

modelCom.indCom.metSps((id + 1): (m + 1)) = modelCom.indCom.metSps(id:m);
modelCom.indCom.metSps(id) = speciesId;
modelCom.indCom.Mcom(modelCom.indCom.Mcom >= id) = modelCom.indCom.Mcom(modelCom.indCom.Mcom >= id) + 1;
modelCom.indCom.Msp(modelCom.indCom.Msp >= id) = modelCom.indCom.Msp(modelCom.indCom.Msp >= id) + 1;
%move the position of met-related fields
field = fieldnames(modelCom);
for j = 1:numel(field)
    if length(field{j}) >= 3 
        if strcmp(field{j}(1:3), 'met') && ~strcmp(field{j}, 'mets') ...
                && ~strcmp(field{j}, 'metNames') && ~strcmp(field{j}, 'metFormulas')
            if size(modelCom.(field{j}), 1) == m
                modelCom.(field{j})((id + 1): (m + 1),:) = modelCom.(field{j})(id:m,:);
                if iscell(modelCom.(field{j})(1,:))
                    modelCom.(field{j})(id,:) = repmat({''}, 1, size(modelCom.(field{j}),2));
                else
                    modelCom.(field{j})(id,:) = 0;
                end
            end
            if size(modelCom.(field{j}), 2) == m
                modelCom.(field{j})(:,(id + 1): (m + 1)) = modelCom.(field{j})(:,id:m);
                if iscell(modelCom.(field{j})(1))
                    modelCom.(field{j})(:,id) = repmat({''}, size(modelCom.(field{j}),1), 1);
                else
                    modelCom.(field{j})(:,id) = 0;
                end
            end
        end
    end
end

end



                
            


