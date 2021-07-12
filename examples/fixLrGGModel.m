%Ethan Green
%February 14th, 2020
%Attempting to fix the LrGG Model

%% Loading in relevant data
clear
load BiGGmetData_01_04_2017.mat
load LrGG_Model.mat

%% General model statistics
%Number of reactions
disp("The LrGG Model has "+length(model.rxns)+" reactions,")
%Number of metabolites
disp(length(model.mets)+" metabolites,")
%Number of genes
disp(length(model.genes)+" genes,")
%Number of transport reactions
rxnEx = sum(model.S ~= 0, 1) <= 1;
disp(sum(rxnEx == 1)+" transport reactions,")
%Number of extracellular metabolites
%metEx =

%%%%%%%%% may need some other processing of the met IDs to match those in the BiGGmetData
mets = model.mets;

%% Fixing metCharges, metFormulas, and mets
[yn, id] = ismember(mets, BiGGmetData(:, 2));
model.metCharges(yn, 1) = str2double(BiGGmetData(id(yn), 7));
model.metFormulas(yn, 1) = BiGGmetData(id(yn), 6);
model.mets = regexprep(model.mets, '\[C_(\w)\]','\[$1\]');
disp("The model's Charges, metabolite names and formulas were standardized using the BiGG Database.")

%% Intervention - adding biomass metabolite
model = addMetabolite(model,'biomass[e]', 'metName','Biomass');
model = addReaction(model,'EX_biomass(e)','reactionFormula','biomass[e] ->','printLevel',0);
% add the biomass metabolite into the biomass reaction
m = findMetIDs(model, 'biomass[e]');
r = findRxnIDs(model, 'Growth');
model.S(m,r) = 1;
disp('A biomass metabolite was added to track growth.')

%% Classifying R and X group containing formulae as unknowns
%getElementalComposition creates to matrices, elements lists all of the
%elements in the metabolite(s) examined and metEle lists the number of
%elements seen in the metabaolite(s) referencing the elements vector
[metEle, elements] = getElementalComposition(model.metFormulas);
%Using pulled data to identify any problematic 'R' or 'X'groups in a
%compound
genEle = ismember(elements, {'R', 'X'});
%Removing all metFormulas identified to have R
l = length(model.metFormulas(any(metEle(:, genEle), 2)));
disp(l+" metabolites were identified with R or X in their formula.")
model.metFormulas(any(metEle(:, genEle), 2)) = {''};

%% Balance check before autocorrect
[model2, metFormulae, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'fillMets', {'HCharge1', 'H2O'});
rxnEx = sum(model.S ~= 0, 1) <= 1;
rxnActive = model.lb ~= 0 | model.ub ~= 0;
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & ~rxnEx & rxnActive');
% check unbalanced reactions
disp(length(rxnImbal)+" reactions were found to be imbalanced.")

%% Filling metabolites
% automatically fill missing protons and H2O
hc = findMetIDs(model, 'h[c]');
h2o = findMetIDs(model, 'h2o[c]');
[metEleFill, elementsFill] = getElementalComposition({'HCharge1', 'H2O'}, elements);
metComp = regexp(model.mets, '\[\w\]$', 'match', 'once');
acceptAutoCorrect = false(numel(model.rxns), 1);
for j = 1:numel(model.rxns)
    % if automatic filling suggested
    if any(S_fill(:, j))
        % check whether it perfectly resolves the imbalance
        if all(abs(metEleFill' * S_fill(:, j) + metEle' * model.S(:, j)) < 1e-5)
            if all(strcmp(metComp(model.S(:, j) ~= 0), '[c]'))  % this is currently trivial since the model is not compartmentalized
                % accept
                acceptAutoCorrect(j) = true;
                model.S([hc; h2o], j) = model.S([hc; h2o], j) + S_fill(:, j);
            end
        end
    end
end

% %% Balance check after autocorrect
% [model2, metFormulae, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'fillMets', {'HCharge1', 'H2O'});
% rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & ~rxnEx & rxnActive');
% % check unbalanced reactions
% disp(length(rxnImbal)+" reactions were found to be imbalanced after automatic intervention.")

%%%%%%%  Add manual curation here to fix unbalanced reactions
% %% Manually Solving Imbalances
% model.S(20,207) = -1; %Reaction 207
% model.metCharges(538) = -1; %Reactions 258, 991
% model.S(20,361) = -1; %Reaction 361
%     model.S(437,361) = 0;
%     model.S(434,361) = 1;
% model.metFormulas(459) = {'C20H38O1Conserve_j'}; %Reaction 436
% model.metCharges(678) = -1; %Reaction 437
%     model.metCharges(679) = -1;
%     model.metCharges(680) = -1;
%     model.metCharges(681) = -1;
% model.metCharges(764) = -1; %Reaction 549
%     model.metCharges(1018) = -1;
% model.metCharges(676) = -1; %Reactions 987, 1033
%     model.metCharges(949) = 0;
%     model.S(20,1033) = -14;
% model.S(20,1216) = 2; %Reaction 1216

%% Balance check manual correction
[model2, metFormulae, elements, metEle, rxnBal, S_fill, solInfo] = computeMetFormulae(model, 'fillMets', {'HCharge1', 'H2O'});
rxnImbal = model.rxns(any(abs(rxnBal) > 1e-4, 1) & ~rxnEx & rxnActive');
% check unbalanced reactions
if length(rxnImbal) > 0
  % This will print the first unbalanced reaction if there is any
    printImbalance(model, rxnImbal(1), 0, rxnBal, elements, metEle)
    disp(length(rxnImbal)+" reactions were found to be imbalanced after manual intervention.")
else
    disp("The model was balanced after manual intervention and is saved as model2.")
end
