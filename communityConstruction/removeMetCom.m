function [modelCom] = removeMetCom(modelCom, metID,metFieldKnown)
%[modelCom] = removeMetCom(modelCom, metID, metFieldKnown)
% Remove metabolites from the model.
% metID is the indices or names of the metabolites to be removed
% Update the field mets, metNames, metFormulas, metSps, b, and other fields
% starting with 'met' and have at least one dimension the same as modelCom.mets
% Other met fields to be updated can be supplied using metFieldKnown
if nargin < 3
    metFieldKnown = {};
end
if ~isnumeric(metID)
    metID = findMetIDs(modelCom,metID);
    if ~all(metID)
        warning('Some metabolite IDs are not present in the model.');
        return
    end
end
m = size(modelCom.S,1);
if isfield(modelCom,'indCom')
    indCom = modelCom.indCom;
else
    indCom = infoCom2indCom(modelCom);
end
%specifically remove rows in EXcom and EXsp if community metabolites are
%deleted
[~, mIdCom] = ismember(metID, find(indCom.metSps == 0));
mIdCom = mIdCom(mIdCom ~= 0);
indCom.EXcom(mIdCom,:) = [];
indCom.EXsp(mIdCom,:) = [];
indCom.Msp(mIdCom,:) = [];
indCom.Msp(ismember(indCom.Msp,metID)) = 0;
Msp = cell(size(indCom.Msp));
Msp(indCom.Msp~=0) = modelCom.mets(indCom.Msp(indCom.Msp~=0));

modelCom.S(metID,:) = [];
field = fieldnames(modelCom);
metField = union(field(strncmp('met',field,3)),metFieldKnown);
metField = setdiff(metField,{'b'});
for j = 1:numel(metField)
    if size(modelCom.(metField{j}), 1) == m
        modelCom.(metField{j})(metID,:) = [];
    end
    if size(modelCom.(metField{j}), 2) == m
        modelCom.(metField{j})(:,metID) = [];
    end
end
if (isfield(modelCom,'b'))
    modelCom.b(metID) = [];
else
    modelCom.b = zeros(length(modelCom.mets),1);
end
indCom.metSps(metID) = [];
indCom.Mcom = find(indCom.metSps == 0);
indCom.Msp(indCom.Msp~=0) = findMetIDs(modelCom, Msp(indCom.Msp~=0));
modelCom.indCom = indCom;
modelCom.infoCom = infoCom2indCom(modelCom,indCom,true,modelCom.infoCom.spAbbr,modelCom.infoCom.spName);

                
            


