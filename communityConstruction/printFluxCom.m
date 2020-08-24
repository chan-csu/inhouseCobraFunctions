function [fluxEXsp, fluxUTsp, fluxEXcom, metId, out] = printFluxCom(model, flux, metFlag, SpFlag, minFlux, printFlag, srFlag, biomass)
%[fluxEXsp, fluxEXcom, metId, out] = printFluxCom(model, flux, metFlag, SpFlag, minFlux, printFlag, srFlag, biomass)
%print a flux vector for community model
%
% input:
% model: the community model with the field EXcom and EXsp indicating
%       indicies of community exchange reactions, and spBm indicating indicies
%       of biomass reaction of each species
% flux: the flux vector to be printed
% SpFlag: true to show metabolites with individual exchange fluxes though 
%         community exchange flux is zero.
%         false (default) not to show.
% metFlag: true to print with model.metNames. 
%          false (default) to print with model.mets
% minFlux: only print metabolites with exchange fluxes > minFlux 
%          (default 1e-6)
% printFlag: true to print (default). false not to print.
% srFlag:   true to print the specific exchange rate (flux/biomass).
%           Default to be false
% biomass:  the biomass vector for printing specific rate. If not given,
%           use flux(model.spBm)
if ~exist('SpFlag', 'var')
    SpFlag = false;
elseif isempty(SpFlag)
    SpFlag = false;
end
if ~exist('metFlag', 'var')
    metFlag = false;
elseif isempty(metFlag)
    metFlag = false;
end
if ~exist('printFlag', 'var')
    printFlag = true;
elseif isempty(printFlag)
    printFlag = true;
end
if ~exist('srFlag', 'var')
    srFlag = false;
elseif isempty(srFlag)
    srFlag = false;
end
if ~isfield(model,'indCom')
    if ~isfield(model,'infoCom') || ~isstruct(model.infoCom) || ...
            ~all(isfield(model.infoCom,{'spBm','EXcom','EXsp','spAbbr','rxnSps','metSps'}))
        error('infoCom must be provided for calculating the max. community growth rate.\n');
    end
    %get useful reaction indices
    model.indCom = infoCom2indCom(model);
end
if isfield(model,'infoCom') && isfield(model.infoCom,'spAbbr')
    model.sps = model.infoCom.spAbbr;
else
    model.sps = strcat('Sp',strtrim(cellstr(num2str((1:size(model.indCom.EXsp,2))'))));
end

if srFlag && ~exist('biomass', 'var')
    %warning('Biomass vector must be supplied for printing specific rates.')
    biomass = flux(model.indCom.spBm);
end
%minimum value of flux for display
if ~exist('minFlux', 'var') || isempty(minFlux)
    minFlux = 1e-6;
end
%the model has separate compartments for (i) community uptake and 
%(ii) community export and inter-organism exchange or not.
UEsep = isfield(model,'UTsp');
% if not separate, check also if there is only one exchange reaction for
% each community metabolite, or one for uptake one for export
if ~UEsep
    if size(model.indCom.EXcom, 2) == 2
        UErxnSep = true;
    elseif size(model.indCom.EXcom, 2) == 1
        UErxnSep = false;
    else
        error('indCom.EXcom can have either one or two columns.')
    end
end
%number of organisms    
nSp = numel(model.indCom.spBm);
%community metabolites
comMet = find(model.indCom.metSps == 0);
if metFlag
    %print community metabolite names
    comMetDisp = model.metNames(comMet);
    for j = 1:numel(comMetDisp)
        if iscell(comMetDisp{j})
            %get the first name only if multiple names exist
            comMetDisp{j} = comMetDisp{j}{1};
        end
    end
else
    %print community metabolite ID
    comMetDisp = cellfun(@(x) strrep(model.mets{x}, '[u]', ''), num2cell(comMet), 'UniformOutput', false);
end

fluxEXsp = zeros(numel(comMet),nSp);
% fluxEXcom = zeros(numel(comMet),2);
if UEsep
    fluxUTsp = zeros(numel(comMet),nSp);
else
    fluxUTsp = [];
end
show = false(numel(comMet), 1);
for j = 1:numel(comMet)
    if SpFlag
        show(j) = show(j) | any(abs(flux(model.indCom.EXsp(j,model.indCom.EXsp(j,:) > 0))) > minFlux);
        if UEsep
            show(j) = show(j) | any(abs(flux(model.indCom.EXsp(j,model.indCom.UTsp(j,:) > 0))) > minFlux);
        end
    end
    show(j) = show(j) | any(abs(flux(model.indCom.EXcom(j,:))) > minFlux);
%     fluxEXcom(j, model.rxnSps(model.EXsp(j,model.EXsp(j,:) > 0))) = flux(model.EXsp(j,model.EXsp(j,:) > 0));
end

flux0 = flux(model.indCom.spBm) == 0;
flux2print = num2cell(full(flux(model.indCom.spBm)));
flux2print(flux0) = {' '};

if printFlag
    len = cellfun(@(x) length(x), comMetDisp);
    fp = sprintf(['%%' num2str(max(len(show)) + ceil(log10(numel(comMet))) + 4) 's  ']);
    printStr0 = cell(nSp + 3, 1);
    printStr0{1} = fp;
    if UEsep
        printStr0(2:3) = {'%-10.3g'};
        printStr0(4:end) = {'%-21.3g'};
        printStr = cell(2*nSp + 3, 1);
        printStr{1} = fp;
        printStr(2:3) = {'%-10.3g'};
        printStr(4:2:end) = {'%-+10.3g'};
        printStr(5:2:end) = {'%-+11.3g'};
    else
        printStr0(2:end) = {'%-10.3g'};
        printStr = printStr0;
    end
    printStr0([false; true(2,1); flux0]) = {'%-10s'};
    if UEsep
        fprintf([fp '%-10s%-10s' repmat('%-21s', 1, nSp) '\n'],...
            'Mets', 'Comm. UT', 'Comm. EX', model.sps{:});
    else
        fprintf([fp repmat('%-10s', 1, nSp + 2) '\n'],...
            'Mets', 'Comm. UT', 'Comm. EX', model.sps{:});
    end
    fprintf([strjoin(printStr0(:)','') '\n'], 'Biomass', '--', '--', flux2print{:});
end
if UEsep
    out = cell(sum(show) + 2, 2*nSp + 3);
    out(1,:) = [{'Mets', 'Comm. UT', 'Comm. EX'} reshape([model.sps(:)';repmat({''},1,nSp)],1,nSp*2)];
    out(2,:) = [{'Biomass', '--', '--'} reshape([flux2print(:)';repmat({''},1,nSp)],1,nSp*2)];
else
    out = cell(sum(show) + 2, nSp + 3);
    out(1,:) = [{'Mets', 'Comm. UT', 'Comm. EX'} model.sps(:)'];
    out(2,:) = [{'Biomass', '--', '--'} flux2print(:)'];
end
k = 2;
for j = 1:numel(comMet)
    if show(j)
        %organism exchange flux
        fluxSp = zeros(nSp,1);
        fluxSp(model.indCom.EXsp(j,:) > 0) = flux(model.indCom.EXsp(j,model.indCom.EXsp(j,:) > 0));
        if srFlag
            fluxSp = fluxSp ./ biomass(:);
        end
        fluxEXsp(j, :) = fluxSp;
        %organism uptake flux (if exist)
        if UEsep
            fluxSp2 = zeros(nSp,1);
            fluxSp2(model.UTsp(j,:) > 0) = flux(model.UTsp(j,model.UTsp(j,:) > 0));
            if srFlag
                fluxSp2 = fluxSp2 ./ biomass(:);
            end
            fluxUTsp(j, :) = fluxSp2;
            %arrange the flux to be [comUT | comEX | -spUT + spEX]
            fluxSp = [fluxSp2 fluxSp]';
        end
        if ~UEsep && ~UErxnSep
            flux2print = full(flux(model.indCom.EXcom(j,:)));
            flux2print = [-min([flux2print, 0]); max([flux2print, 0]); fluxSp(:)];
        else
            %formatting
            flux2print = [full(flux(model.indCom.EXcom(j,:))); fluxSp(:)];
        end
        flux0 = abs(flux2print) < minFlux;
        flux2print = num2cell(flux2print);
        flux2print(flux0) = {' '};
        k = k + 1;
        out{k, 1} = ['(' num2str(j) ') ' comMetDisp{j}];
        out(k, 2:end) = flux2print'; 
        if printFlag
            printStrJ = printStr;
            if UEsep
                UTlogic = reshape([true(1,nSp); false(1,nSp)],2*nSp,1);
                printStrJ([false; flux0(1:2); false(2*nSp,1)]) = {'%-10s'};
                printStrJ([false(3,1); flux0(3:end) & UTlogic]) = {'%-10s'};
                printStrJ([false(3,1); flux0(3:end) & ~UTlogic]) = {'%-11s'};
            else
                printStrJ([false; flux0]) = {'%-10s'};
            end
            fprintf([strjoin(printStrJ(:)','') '\n'], out{k,:});
        end
    end
end
fluxEXcom = flux(model.indCom.EXcom);
metId = find(show);

end

