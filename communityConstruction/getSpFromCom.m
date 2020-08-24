function [modelCell, info] = getSpFromCom(modelCom, SpAbbr, exAbbr, rxnIden, metIden,rxnFromMet,metFromRxn)
%[modelCell, info] = getSpFromCom(modelCom, SpAbbr, exAbbr, rxnIden, metIden)
%Get separate species COBRA model from a joint community model
%Input:
%  modelCom:   Community model
%  SpAbbr:     Cell array of organisms' abbreviation used in modelCom.rxns
%              and modelCom.mets as unique identifiers
%  exAbbr:     Identifier for community metabolites (default '[u]')
%  rxnIden:    A regular expression for extracting the reaction names, with 
%              a '%s' denoting the organism's abbreviation in SpAbbr, e.g.,
%              '(_%s)$'. The token in the parenthesis will be removed. 
%              If a reaction is named 'ABC_Ec' in organism Ec,
%              'ABC' would be extracted. (default '(_%s)$')
%  rxnIden:    A regular expression for extracting the metabolite names, with 
%              a '%s' denoting the organism's abbreviation in SpAbbr.
%              (default '(_%s)\]$')
%  rxnFromMet  true to identify rxns from identified mets for each organism
%  metFromRxn  true to identify mets from identified rxns for each organism

if ~exist('SpAbbr', 'var') || isempty(SpAbbr)
    spNameFromModel = true;
else
    spNameFromModel = false;
end
if spNameFromModel
    if isfield(modelCom, 'infoCom') && isfield(modelCom.infoCom,'spAbbr')
        SpAbbr = modelCom.infoCom.spAbbr;
    else
        error('Abbreviation of organisms must be supplied.')
    end
end
if ~exist('exAbbr', 'var') || isempty(exAbbr)
    exAbbr = '[u]';
end
nSp = numel(SpAbbr);
metrxnFromCom = false;
if isfield(modelCom,'indCom') || isfield(modelCom,'infoCom')
    if isfield(modelCom,'indCom')
        indCom = modelCom.indCom;
    else isfield(modelCom,'infoCom');
        indCom = infoCom2indCom(modelCom);
    end
    %mets and rxns for each organism given
    metrxnFromCom = true;
end
if ~exist('metIden','var') || isempty(metIden)
    metIden = '(_%s)\]$';
end
if ~exist('rxnIden','var') || isempty(rxnIden)
    rxnIden = '(_%s)$';
end
if ~exist('rxnFromMet','var') || isempty(rxnFromMet)
    rxnFromMet = false;
end
if ~exist('metFromRxn','var') || isempty(metFromRxn)
    metFromRxn = ~metrxnFromCom;
end
if metFromRxn && rxnFromMet && ~metrxnFromCom
    warning('No information distinguishing different organisms provided. Turn on rxnFromMet.')
    rxnFromMet = false;
end

[metSpAbbr,rxnSpAbbr]  = deal(repmat({''},nSp,1));
for jSp = 1:nSp
    metSpAbbr{jSp} = strrep(metIden,'%s',SpAbbr{jSp});
    rxnSpAbbr{jSp} = strrep(rxnIden,'%s',SpAbbr{jSp});
end

modelCell = num2cell(repmat(modelCom, nSp, 1));
info = struct();

[m, n] = size(modelCom.S);

%species-specific reactions
rxnSp = false(n, nSp);
metSp = false(m, nSp);
if metrxnFromCom
    for j = 1:nSp
        rxnSp(:,j) = indCom.rxnSps == j;
        metSp(:,j) = indCom.metSps == j;
    end
else
    if ~rxnFromMet
        for j = 1:nSp
            rxnSp(:,j) = ~cellfun(@isempty, regexp(modelCom.rxns,rxnSpAbbr{j},'once'));
        end
    end
    if ~metFromRxn
        for j = 1:nSp
            metSp(:,j) = ~cellfun(@isempty, regexp(modelCom.rxns,metSpAbbr{j},'once'));
        end
    end
    if rxnFromMet
        for j = 1:nSp
            rxnSp(:,j) = any(modelCom.S(metSp(:,j), :), 1);
        end
    end
    if metFromRxn
        for j = 1:nSp
            metSp(:,j) = any(modelCom.S(:, rxnSp(:, j)), 2);
        end
    end
end

if any(sum(rxnSp, 2) > 1)
    warning('The provided identifiers are not unique in some reactions. Return.')
    modelCell = {};
    info.rxnNotUnique = modelCom.rxns(sum(rxnSp, 2) > 1);
    return
end
if any(sum(metSp, 2) > 1)
    warning('The provided identifiers are not unique in some metabolites. Return.')
    modelCell = {};
    info.metNotUnique = modelCom.mets(sum(metSp, 2) > 1);
    return
end

%community exchange reactions
comRxnEx = ~any(rxnSp, 2);
%possibly problematic community reactions
comRxnProblem = comRxnEx & sum(modelCom.S ~= 0)' ~= 1;
info.comRxnProblem = find(comRxnProblem);
comRxnEx = comRxnEx & ~comRxnProblem;

%community metabolites found by exchange reactions
comMet1 = any(modelCom.S(:, comRxnEx), 2);
%alternatively found by identifier '[u]'
comMet2 = cellfun(@(x) ~isempty(strfind(x, exAbbr)), modelCom.mets);
if ~isequal(comMet1, comMet2)
    %conflict in metabolites
    info.comMetNoEx = comMet2 & ~comMet1;
    info.comMetNoU = comMet1 & ~comMet2;
end

%assign different fields
fields = fieldnames(modelCom);
%fields = setdiff(fieldnames(modelCom), {'S'});
for j = 1:numel(fields)
    if (size(modelCom.(fields{j}), 1) == m && size(modelCom.(fields{j}), 2) == n)
        %S matrix
        flag = 1;
    elseif (size(modelCom.(fields{j}), 1) == n && size(modelCom.(fields{j}), 2) == m)
        %rxn x met
        flag = 2;
    elseif size(modelCom.(fields{j}), 1) == n
        %#rxn rows
        flag = 3;
    elseif size(modelCom.(fields{j}), 2) == n
        %#rxn coloumns
        flag = 4;
    elseif size(modelCom.(fields{j}), 1) == m
        %#met rows
        flag = 5;
    elseif size(modelCom.(fields{j}), 2) == m
        %#met coloumns
        flag = 6;
    else
        flag = 0;
    end
    if strcmp(fields{j},'genes')
        %special care for genes
        flag = -1;
    end
    for k = 1:nSp
        if strcmp(fields{j},'rxnGeneMat')
            %special care for rxnGeneMat
            geneK = any(modelCom.rxnGeneMat(rxnSp(:,k),:),1);
            modelCell{k}.rxnGeneMat = modelCom.rxnGeneMat(rxnSp(:,k),geneK);
            modelCell{k}.genes = modelCom.genes(geneK);
        elseif flag == 1
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(metSp(:, k), rxnSp(:,k));
        elseif flag == 2
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(rxnSp(:, k), metSp(:,k));
        elseif flag == 3
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(rxnSp(:, k), :);
        elseif flag == 4
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(:, rxnSp(:, k));
        elseif flag == 5
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(metSp(:, k), :);
        elseif flag == 6
            modelCell{k}.(fields{j}) = modelCom.(fields{j})(:, metSp(:, k));
        elseif flag == 0
            %if none of the above, just copy the whole from the joint model
            modelCell{k}.(fields{j}) = modelCom.(fields{j});
        end
        if strcmp(fields{j}, 'rxns') || strcmp(fields{j}, 'mets')
            %remove the species identifier for strings
            nameEmpty = cellfun(@isempty,modelCell{k}.(fields{j}));
            modelCell{k}.(fields{j})(nameEmpty) = repmat({''}, sum(nameEmpty), 1);
            if strcmp(fields{j}, 'rxns')
                str = rxnSpAbbr{k};
            else
                str = metSpAbbr{k};
            end
            position = regexp(modelCell{k}.(fields{j}),str,'tokenExtents');
            %position = ~cellfun(@isempty,position);
            %first check if any tokens are in metIden or rxnIden
            token = false;
            for jF = 1:numel(modelCell{k}.(fields{j}))
                if ~isempty(position{jF})
                    for jC = 1:numel(position{jF})
                        if ~isempty(position{jF}{jC})
                            token = true;
                            break
                        end
                    end
                end
            end
            if token
                %there are tokens
                for jF = 1:numel(modelCell{k}.(fields{j}))
                    if ~isempty(position{jF})
                        str2 = true(1,length(modelCell{k}.(fields{j}){jF}));
                        for jC = 1:numel(position{jF})
                            if ~isempty(position{jF}{jC})
                                str2(position{jF}{jC}(1):position{jF}{jC}(2)) = false;
                            end
                        end
                        modelCell{k}.(fields{j}){jF} = modelCell{k}.(fields{j}){jF}(str2);
                    end
                end
            else
                %no tokens. Directly replace everything
                modelCell{k}.(fields{j}) = regexprep(modelCell{k}.(fields{j}),str,'');
            end
        end
    end
    if flag == 1
        info.(fields{j}) = modelCom.(fields{j})(comMet2, comRxnEx);
    elseif flag == 2
        info.(fields{j}) = modelCom.(fields{j})(comRxnEx, comMet2);
    elseif flag == 3
        info.(fields{j}) = modelCom.(fields{j})(comRxnEx, :);
    elseif flag == 4
        info.(fields{j}) = modelCom.(fields{j})(:, comRxnEx);
    elseif flag == 5
        info.(fields{j}) = modelCom.(fields{j})(comMet2, :);
    elseif flag == 6
        info.(fields{j}) = modelCom.(fields{j})(:, comMet2);        
    end
end

%get link between community metabolites and species-specific metabolites
spExProblem = [];
for j = 1:nSp
    modelCell{j}.metComs = cell(sum(metSp(:,j)), 1);
    %species-community exchange reactions
    exJ = find(sum(modelCom.S(comMet2,:) ~= 0)' > 0 & rxnSp(:,j));
    for k = 1:numel(exJ)
        id = find(modelCom.S(:, exJ(k)));
        if numel(id) ~= 2
            spExProblem = [spExProblem; exJ(k)];
        else
            if metSp(id(2),j)
                id = id([2 1]);
            end
            if metSp(id(1),j) && ~metSp(id(2),j)
                modelCell{j}.metComs{sum(metSp(1:id(1),j))} = modelCom.mets{id(2)};
            else
                spExProblem = [spExProblem; exJ(k)];
            end
        end
    end
end
info.spExProblem = spExProblem;