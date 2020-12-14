%% data and parameters needed
% directory of the models
modelDir = 'modelsChan2017';
% model files to load
modelsToLoad = {'iAH991_norm.mat'; ...  Bacteroides thetaiotamicron VPI-5482. Heinken et al, 2013. Updated in Heinken & Thiele, 2015a.	
    'iEre400_norm.mat'; ...  Eubacterium rectale ATCC33656. Shoaie et al, 2013
    'iFpraus_v1_norm.mat'; ...  Faecalibacterium prausnitzii A2-165. Heinken et al, 2014
    'Ef_V583_norm.mat'; ...  Enterococcus faecalis V583.	Veith et al, 2015
    'iLca12A_norm.mat'; ...  Lactobacillus casei ATCC 334	Vinay-Lara et al, 2014
    'iMP429_norm.mat'; ...  Streptococcus thermophilus LMG18311	Pastink et al, 2009, Updated in Heinken & Thiele, 2015a.
    'iBif452_norm.mat'; ...  Bifidobacterium adolescentis L2-32	El-Semman et al, 2014
    'iJO1366_update.mat'; ...  Escherichia coli K-12 MG1655	Orth et al, 2011
    'iYL1228_norm.mat'};  % Klebsiella pneumonia MGH 78578	Liao et al, 2011. Updated in Heinken & Thiele, 2015a.

% Info that helps define the index matrices for the community model when
% provided. See help createCommModel
options = struct();
% biomass reaction name
options.spBm = {'Biomass_BT_v2_norm','r55_norm','Biomass_FP_norm','BIOMASS_norm', ...
    'bio02514_norm','biomass_norm_STR','Biomass_norm','Ec_biomass_iJO1366_WT_53p95M','Biomass_norm'};
% species abbreviation
options.spAbbr = {'Bt', 'Er', 'Fp', 'Ef', 'Lc', 'St', 'Ba', 'Ec', 'Kp'};
% ATPM reaction
options.spATPM = {'ATPM', 'r54', 'ATPM', 'ATPM', 'ATPM', 'ATPM', 'ATPM', 'ATPM', 'ATPM'};
% single system exchange reaction for each community metabolite
options.sepUtEx = false;
% detect extracellular metabolites using [e]
options.metExId = '[e]';

% ModelSeed data for mapping to BiGG IDs
seedMet = load('seedMet201907.mat');
seedMet = seedMet.seedMet;
seedRxn = load('seedRxn201907.mat');
seedRxn = seedRxn.seedRxn;

%% Load the model and unify the IDs if needed
models = cell(size(modelsToLoad));

for j = 1:numel(modelsToLoad)
    % Can't read those model using readCbModel since some formatting of the
    % model structures are already outdated. Need some cosmetic fixing
    %     models{j} = readCbModel([modelDir filesep modelsToLoad{j}]);
    % Use load model directly instead
    model = load([modelDir filesep modelsToLoad{j}]);
    models{j} = model.model;
end
clear model

for j = [3, 4, 6, 9]
    models{j}.metBiGG((end + 1):numel(models{j}.mets)) = {''};
end

for j = 1:numel(models)
    if sum(strncmp(models{j}.mets, 'cpd', 3)) / numel(models{j}.mets) > 0.5
        % more than half starting with CPD, assuming Seed ID
        % Conversion to BiGG ID
        seedID = regexp(models{j}.mets, 'cpd\d{5}', 'match', 'once');
        seedID(cellfun(@isempty, seedID)) = {''};
        [yn, id] = ismember(seedID, seedMet(:, 1));
        for k = 1:numel(models{j}.mets)
            if yn(k)
                biggID = regexp(seedMet{id(k), 19}, 'BiGG\:([^;]+);?', 'tokens', 'once');
                if ~isempty(biggID)
                    biggID = biggID{1};
                
                    if any(strfind(biggID, '|'))
                        biggID = strsplit(biggID, '|');
                        
                        % Can add the following manual selection here if there are multiple
                        % BiGG IDs:
                        % %      fprintf('%s\t%s\t%s\n', models{j}.mets{k}, models{j}.metNames{k}, models{j}.metFormulas{k})
                        % %      for p = 1:numel(biggID)
                        % %          fprintf('%d\t%s\n', p, biggID{p})
                        % %      end
                        % %      choice = input('Which one? Input an integer: ', 's');
                        % %      while str2double(choice) < 1 || str2double(choice) > numel(biggID) || floor(str2double(choice)) ~= str2double(choice)
                        % %          choice = input('Which one? Input an integer: ', 's');
                        % %      end
                        % %      biggID = biggID{str2double(choice)};
                        biggID = biggID{1};
                    end
                    
                    models{j}.mets{k} = strrep(models{j}.mets{k}, seedID{k}, biggID);
                end
            end
        end
        
        seedID = regexp(models{j}.rxns, 'rxn\d{5}', 'match', 'once');
        seedID(cellfun(@isempty, seedID)) = {''};
        [yn, id] = ismember(seedID, seedRxn(:, 1));
        for k = 1:numel(models{j}.rxns)
            if yn(k)
                biggID = regexp(seedRxn{id(k), 13}, 'BiGG\:([^;]+);?', 'tokens', 'once');
                if ~isempty(biggID)
                    biggID = biggID{1};
                
                    if any(strfind(biggID, '|'))
                        biggID = strsplit(biggID, '|');
                        % Can add the following manual selection here if there are multiple
                        % BiGG IDs:
                        % %      fprintf('%s\t%s\t%s\n', models{j}.mets{k}, models{j}.metNames{k}, models{j}.metFormulas{k})
                        % %      for p = 1:numel(biggID)
                        % %          fprintf('%d\t%s\n', p, biggID{p})
                        % %      end
                        % %      choice = input('Which one? Input an integer: ', 's');
                        % %      while str2double(choice) < 1 || str2double(choice) > numel(biggID) || floor(str2double(choice)) ~= str2double(choice)
                        % %          choice = input('Which one? Input an integer: ', 's');
                        % %      end
                        % %      biggID = biggID{str2double(choice)};
                        biggID = biggID{1};
                    end
                
                    models{j}.rxns{k} = strrep(models{j}.rxns{k}, seedID{k}, biggID);
                end
            end
        end
    elseif isfield(models{j}, 'metBiGG')
        if ~all(cellfun(@(x) strcmp(x(end), ']'), models{j}.mets))
            error('Check compartment notations')
        end
        
        for k = 1:numel(models{j}.mets)
            if ~isempty(models{j}.metBiGG{k})
                models{j}.mets{k} = regexprep(models{j}.mets{k}, '^[^\[]+(\[[^\]]*\])$', [models{j}.metBiGG{k}, '$1']);
            end
        end
    end
    fn = fieldnames(models{j});
    for k = 1:numel(fn)
        % match the dimensions of some fields in the model
        if strncmp(fn{k}, 'met', 3) && size(models{j}.(fn{k}), 1) < numel(models{j}.mets)
            if isnumeric(models{j}.(fn{k}))
                models{j}.(fn{k})((end + 1):numel(models{j}.mets)) = NaN;
            elseif iscell(models{j}.(fn{k}))
                models{j}.(fn{k})((end + 1):numel(models{j}.mets)) = {''};
            end
        end
        if iscell(models{j}.(fn{k}))
            models{j}.(fn{k})(cellfun(@isempty, models{j}.(fn{k}))) = {''};
        end
    end
    models{j}.csense = repmat('E', numel(models{j}.mets), 1);
end

%% compile community model

modelCom = createCommModel(models, options);
% check out modelCom.infoCom and modelCom.indCom. They are useful in
% manipulating the model