% load the Clostridium ljungdahlii model
model = readCbModel('12934_2013_907_MOESM6_ESM.mat');

% we would like to check if there is any internal ATP cycle such that
% without consuming any substrate, the cell can generate ATP.
% For example, if there are reactions:
%
% R1: a[c] + adp[c] -> b[c] + atp[c]      (hypothetical reaction)
% R2: b[c] + pi[c] + h[c] -> a[c] + h2o[c]  (hypothetical reaction)
% ATPM: atp[c] + h2o[c] -> adp[c] + pi[c] + h[c]  (ATP maintenance reaction/ATP hydrolysis)
% (pi[c] is cytoplasmic phosphate and h[c] is cytoplasmic proton, or hydrogen ion)
%
% Then [v_R1, v_R2, v_ATPM] = [1 1 1] satisfies Sv = 0 (mass balance),
% which means some intracellular reactions alone can give rise to ATP
% hydrolysis. So this is problematic.
% To see if this exists, we first shut down all uptake reactions and then
% maximize ATPM. If you get nonzero flux, there are problems.

% first make a copy of the model and change the model objective to ATPM flux
model2 = changeObjective(model, 'ATPM');

% get all exchange reactions (recall: all columns with only one nonzero stoichiometry)
rxnEx = sum(model.S ~= 0, 1) == 1;
% change all LBs to 0
model2.lb(rxnEx) = 0;

% run FBA
% - 'max' means maximizing the objective function .c (default is max already)
% - 'one' means finding a flux vector that not only maximizes biomass, but
%   also minimizes the total fluxes. Easier to look up if there are ATP cycles
s = optimizeCbModel(model2, 'max', 'one');
% infeasible in this case (s.stat = 0, 'stat' for solution status)
% ('assert' returns error if the statement is not true. For checking purpose)
assert(s.stat == 0)

% It is because the lower bound for ATPM is positve (0.45):
surfNet(model2, 'ATPM')
% Therefore when the model cannot have any uptake, it cannot generate any
% ATP so it cannot satisfy the positive ATPM requirement.

% However, if we don't impose a +ve lower bound for ATPM:
model2 = changeRxnBounds(model2, 'ATPM', 0, 'l');
surfNet(model2, 'ATPM')

% If we run FBA again:
s = optimizeCbModel(model2, 'max', 'one');
% now optimal solution can be found (s.stat = 1) but the objective function
% value is 0 (s.f = 0), which implies that the model cannot make any ATP in
% the absense of any uptake
assert(s.stat == 1 & s.f == 0)

%% Let's make an ATP cycle!!!

% PIabc is the ABC transporter for phosphate, which consumes ATP to
% transport, a type of active transport.
surfNet(model2, 'PIabc')

% Now if we make it reversible:
model2 = changeRxnBounds(model2, 'PIabc', -1000, 'l');
surfNet(model2, 'PIabc')

% now the lower bound is -1000. Run FBA again:
s2 = optimizeCbModel(model2, 'max', 'one');
% We get positive objective function value, which means ATP can be produced
% while there is still no uptake.
assert(s2.stat == 1 & s2.f > 0)

% We want to check what the problem is. Look at the active reactions
surfNet(model2, model2.rxns(s2.x ~= 0), [], s2.x, 1, 0) 
% s2.x ~= 0 is the logical index vector for reactions with nonzero fluxes
% model2.rxns(s2.x ~= 0) gets the names for those reactions
%
% the last three inputs 's2.x, 1, 0)' represent:
% - the flux vector (s2.x)
% - show nonzero fluxes only (1)
% - not showing details of each metabolite in a reaction (0)
%
% You should see four reactions shown: ATPM, ATPS4r, PIabc, PIt2r with
% fluxes 750, -250, -1000, 1000 respective.
% So basically:
% - PIabc keeps pumping phosphate out of the cell but at the same time gains 1000 ATP
% - PIt2r recyle the 1000 phosphate in the cost of 1000 extracellular proton
% - ATPS4r uses 250 ATP and pumps out 1000 extracellular protons 
% - There are 750 ATP left for ATP hydrolysis (ATPM)
%
% In this case, you should be able to figure out that PIabc being
% reversible does not make sense, so change it to irreversible. Then there
% is no ATP cycle.