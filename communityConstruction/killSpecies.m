function [modelK] = killSpecies(model, species)
%[modelK] = killSpecies(model, species)
%return a model shutting off all reactions of the given species
%(abbreviations or Ids)
if ischar(species)
    species = {species};
end

if iscell(species)
    del = false(numel(species), 1);
    species2 = zeros(numel(species),1);
    for j = 1:numel(species)
        f = find(strcmp(species{j}, model.sps));
        if ~isempty(f)
            species2(j) = f;
        elseif ~strcmp(species{j}, 'com') 
            del(j) = true;
        end
        %else if species{j} = 'com', not allowing community exchange, for checking loops
    end
    species2(del) = [];
elseif (islogical(species) || all(species == 0 | species == 1)) && numel(species) == numel(model.indCom.spBm)
    species2 = find(species);
elseif isnumeric(species)
    species2 = species(species <= numel(model.indCom.spBm));
else
    warning('Incorrect input species.')
    modelK = [];
    return
end


modelK = model;
for j = 1:numel(species2)
    modelK.lb(modelK.indCom.rxnSps == species2(j)) = 0;
    modelK.ub(modelK.indCom.rxnSps == species2(j)) = 0;
end
end

