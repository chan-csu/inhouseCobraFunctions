function [spName,spId] = activeSpecies(modelCom, flux)
%spInd = activeSpecies(modelCom, flux)
%
%Given a flux vector 'flux' and a community model 'modelCom' with the field
%'rxnSps'. Return the IDs of species that are active in the flux
%distribution.
%
%If a flux distribution is not provided, return the names and IDs of the
%species with any non-zero bounds

if nargin == 2
    spId = setdiff(unique(modelCom.indCom.rxnSps(abs(flux) > 1e-7)), 0);
elseif nargin == 1
    spId = false(numel(modelCom.infoCom.spAbbr),1);
    for j = 1:numel(modelCom.infoCom.spAbbr)
        spId(j) = any(modelCom.ub(modelCom.indCom.rxnSps == j)) | any(modelCom.lb(modelCom.indCom.rxnSps == j));
    end
    spId = find(spId);
end
spName = modelCom.infoCom.spAbbr(spId);