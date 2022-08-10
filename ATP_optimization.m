clear 
pathModel = 'YOUR_PATH';
filename = 'C08_seed.mat';
load([pathModel, filename])
model = C08_seed;
clear C08_seed

model.osenseStr = char('max');

metaboliteMature= {'NADPH[c]'};
model = addReaction(model,'Sink_NADPH[c]', 'reactionFormula', [metaboliteMature{1,1}, ' <=>']);


%% 28 DAF

model28 = model;

% ATP Maintenance
ATPMaintenance = find(ismember(model28.rxns,'RXN-11109_c'));

Biomass = find(ismember(model28.rxns,'BIOMASS-RXN_C08_28DAF'));
GrowthPerHour = 0.306275525/24;
model28.lb(Biomass) = GrowthPerHour;
model28.ub(Biomass) = GrowthPerHour;

% Nutrients according 28DAF:
rxnID = findRxnIDs(model28,'Exchange_Sucrose');
model28.lb(rxnID,1) = -41;
rxnID = findRxnIDs(model28,'Exchange_GLC');
model28.lb(rxnID,1) = -5;
rxnID = findRxnIDs(model28,'Exchange_GLN');
model28.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model28,'Exchange_ASN');
model28.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model28,'Exchange_O2');
model28.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model28,'Exchange_Photon');
model28.lb(rxnID,1) = -10;


% only one active Biomass rxn
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_40DAF');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_60DAF');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_MATURE');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_28DAF_ISOFLAVONES');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_40DAF_ISOFLAVONES');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_60DAF_ISOFLAVONES');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model28,'BIOMASS-RXN_C08_MATURE_ISOFLAVONES');
model28.lb(rxnID,1) = 0;
model28.ub(rxnID,1) = 0;


% OBJ Function
model28 = changeObjective(model28,'RXN-11109_c',1);


% Add ratios
% here define your ratio
model28_Ratio = addRatioReaction(model28,{'RXN-11109_c' 'Sink_NADPH[c]'}, [1 15]);


ATP_FBAsolution28 = optimizeCbModel(model28_Ratio);


fprintf('\n28DAF conditions ATP max using ATP:NADPH ratio = 1:15 : %.4f\n', ATP_FBAsolution28.f);
Biomass28 = findRxnIDs(model28_Ratio,'BIOMASS-RXN_C08_28DAF');
fprintf('Biomass without ISO flux: %.4f\n', ATP_FBAsolution28.x(Biomass28));
Biomass28ISO = findRxnIDs(model28_Ratio,'BIOMASS-RXN_C08_28DAF_ISOFLAVONES');
fprintf('Biomass with ISO flux: %.4f\n', ATP_FBAsolution28.x(Biomass28ISO));
Sink_NADPH = findRxnIDs(model28_Ratio,'Sink_NADPH[c]');
fprintf('Sink_NADPH[c] flux: %.4f\n', ATP_FBAsolution28.x(Sink_NADPH));
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', ATP_FBAsolution28.x(ISO28));
ISO40 = findRxnIDs(model28_Ratio,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', ATP_FBAsolution28.x(ISO40));
ISO60 = findRxnIDs(model28_Ratio,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', ATP_FBAsolution28.x(ISO60));
ISOM = findRxnIDs(model28_Ratio,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', ATP_FBAsolution28.x(ISOM));



%% 40 DAF
model40 = model;

% ATP Maintenance
ATPMaintenance = find(ismember(model40.rxns,'RXN-11109_c'));


% Add biomass  constraint
Biomass = find(ismember(model40.rxns,'BIOMASS-RXN_C08_40DAF'));
GrowthPerHour = 0.471988925/24;
model40.lb(Biomass) = GrowthPerHour;
model40.ub(Biomass) = GrowthPerHour;

% Nutrients according 40DAF:
rxnID = findRxnIDs(model40,'Exchange_Sucrose');
model40.lb(rxnID,1) = -41;
rxnID = findRxnIDs(model40,'Exchange_GLC');
model40.lb(rxnID,1) = -5;
rxnID = findRxnIDs(model40,'Exchange_GLN');
model40.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model40,'Exchange_ASN');
model40.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model40,'Exchange_O2');
model40.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model40,'Exchange_Photon');
model40.lb(rxnID,1) = -10;

% only one active Biomass rxn
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_28DAF');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
% rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_40DAF');
% model40.lb(rxnID,1) = 0;
% model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_60DAF');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_MATURE');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_28DAF_ISOFLAVONES');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_40DAF_ISOFLAVONES');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_60DAF_ISOFLAVONES');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model40,'BIOMASS-RXN_C08_MATURE_ISOFLAVONES');
model40.lb(rxnID,1) = 0;
model40.ub(rxnID,1) = 0;


% OBJ Function
model40 = changeObjective(model40,'RXN-11109_c',1);


% Add ratios
model40_Ratio = addRatioReaction(model40,{'RXN-11109_c' 'Sink_NADPH[c]'}, [1 15]);

ATP_FBAsolution40 = optimizeCbModel(model40_Ratio);


fprintf('\n40DAF conditions ATP max using ATP:NADPH ratio = 1:15 : %.4f\n', ATP_FBAsolution40.f);
Biomass40 = findRxnIDs(model40_Ratio,'BIOMASS-RXN_C08_40DAF');
fprintf('Biomass without ISO flux: %.4f\n', ATP_FBAsolution40.x(Biomass40));
Biomass40ISO = findRxnIDs(model40_Ratio,'BIOMASS-RXN_C08_40DAF_ISOFLAVONES');
fprintf('Biomass with ISO flux: %.4f\n', ATP_FBAsolution40.x(Biomass40ISO));
Sink_NADPH = findRxnIDs(model40_Ratio,'Sink_NADPH[c]');
fprintf('Sink_NADPH[c] flux: %.4f\n', ATP_FBAsolution40.x(Sink_NADPH));
ISO28 = findRxnIDs(model40_Ratio,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', ATP_FBAsolution40.x(ISO28));
ISO40 = findRxnIDs(model40_Ratio,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', ATP_FBAsolution40.x(ISO40));
ISO60 = findRxnIDs(model40_Ratio,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', ATP_FBAsolution40.x(ISO60));
ISOM = findRxnIDs(model40_Ratio,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', ATP_FBAsolution40.x(ISOM));

%% 60 DAF
model60 = model;

% ATP Maintenance
ATPMaintenance = find(ismember(model60.rxns,'RXN-11109_c'));


% Add biomass  constraint
Biomass = find(ismember(model60.rxns,'BIOMASS-RXN_C08_60DAF'));
GrowthPerHour = 0.258666711/24;
model60.lb(Biomass) = GrowthPerHour;
model60.ub(Biomass) = GrowthPerHour;

% Nutrients according 60DAF:
rxnID = findRxnIDs(model60,'Exchange_Sucrose');
model60.lb(rxnID,1) = -31;
rxnID = findRxnIDs(model60,'Exchange_GLC');
model60.lb(rxnID,1) = -4;
rxnID = findRxnIDs(model60,'Exchange_GLN');
model60.lb(rxnID,1) = -14;
rxnID = findRxnIDs(model60,'Exchange_ASN');
model60.lb(rxnID,1) = -14;
rxnID = findRxnIDs(model60,'Exchange_O2');
model60.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model60,'Exchange_Photon');
model60.lb(rxnID,1) = -1;

% only one active Biomass rxn
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_28DAF');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_40DAF');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
% rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_60DAF');
% model60.lb(rxnID,1) = 0;
% model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_MATURE');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_28DAF_ISOFLAVONES');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_40DAF_ISOFLAVONES');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_60DAF_ISOFLAVONES');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;
rxnID = findRxnIDs(model60,'BIOMASS-RXN_C08_MATURE_ISOFLAVONES');
model60.lb(rxnID,1) = 0;
model60.ub(rxnID,1) = 0;


% OBJ Function
model60 = changeObjective(model60,'RXN-11109_c',1);


% Add ratios
model60_Ratio = addRatioReaction(model60,{'RXN-11109_c' 'Sink_NADPH[c]'}, [1 3]);

ATP_FBAsolution60 = optimizeCbModel(model60_Ratio);


fprintf('\n60DAF conditions ATP max using ATP:NADPH ratio = 1:3 : %.4f\n', ATP_FBAsolution60.f);
Biomass60 = findRxnIDs(model60_Ratio,'BIOMASS-RXN_C08_60DAF');
fprintf('Biomass without ISO flux: %.4f\n', ATP_FBAsolution60.x(Biomass60));
Biomass60ISO = findRxnIDs(model60_Ratio,'BIOMASS-RXN_C08_60DAF_ISOFLAVONES');
fprintf('Biomass with ISO flux: %.4f\n', ATP_FBAsolution60.x(Biomass60ISO));
Sink_NADPH = findRxnIDs(model60_Ratio,'Sink_NADPH[c]');
fprintf('Sink_NADPH[c] flux: %.4f\n', ATP_FBAsolution60.x(Sink_NADPH));
ISO28 = findRxnIDs(model60_Ratio,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', ATP_FBAsolution60.x(ISO28));
ISO40 = findRxnIDs(model60_Ratio,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', ATP_FBAsolution60.x(ISO40));
ISO60 = findRxnIDs(model60_Ratio,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', ATP_FBAsolution60.x(ISO60));
ISOM = findRxnIDs(model60_Ratio,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', ATP_FBAsolution60.x(ISOM));

%% Mature

modelM = model;

% ATP Maintenance 
ATPMaintenance = find(ismember(modelM.rxns,'RXN-11109_c'));


% Add biomass  constraint
Biomass = find(ismember(modelM.rxns,'BIOMASS-RXN_C08_MATURE'));
GrowthPerHour = 0.223810871/24;
modelM.lb(Biomass) = GrowthPerHour;
modelM.ub(Biomass) = GrowthPerHour;

% Nutrients according Mature:
rxnID = findRxnIDs(modelM,'Exchange_Sucrose');
modelM.lb(rxnID,1) = -15;
rxnID = findRxnIDs(modelM,'Exchange_GLC');
modelM.lb(rxnID,1) = -2;
rxnID = findRxnIDs(modelM,'Exchange_GLN');
modelM.lb(rxnID,1) = -17;
rxnID = findRxnIDs(modelM,'Exchange_ASN');
modelM.lb(rxnID,1) = -17;
rxnID = findRxnIDs(modelM,'Exchange_O2');
modelM.lb(rxnID,1) = -1;
rxnID = findRxnIDs(modelM,'Exchange_Photon');
modelM.lb(rxnID,1) = -1;

% only one active Biomass rxn
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_28DAF');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_40DAF');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_60DAF');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
% rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_MATURE');
% modelM.lb(rxnID,1) = 0;
% modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_28DAF_ISOFLAVONES');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_40DAF_ISOFLAVONES');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_60DAF_ISOFLAVONES');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;
rxnID = findRxnIDs(modelM,'BIOMASS-RXN_C08_MATURE_ISOFLAVONES');
modelM.lb(rxnID,1) = 0;
modelM.ub(rxnID,1) = 0;


% OBJ Function
modelM = changeObjective(modelM,'RXN-11109_c',1);


% Add ratios
modelM_Ratio = addRatioReaction(modelM,{'RXN-11109_c' 'Sink_NADPH[c]'}, [1 3]);

ATP_FBAsolutionM = optimizeCbModel(modelM_Ratio);


fprintf('\nMature conditions ATP max using ATP:NADPH ratio = 1:3 : %.4f\n', ATP_FBAsolutionM.f);
BiomassM = findRxnIDs(modelM_Ratio,'BIOMASS-RXN_C08_MATURE');
fprintf('Biomass without ISO flux: %.4f\n', ATP_FBAsolutionM.x(BiomassM));
BiomassMISO = findRxnIDs(modelM_Ratio,'BIOMASS-RXN_C08_MATURE_ISOFLAVONES');
fprintf('Biomass with ISO flux: %.4f\n', ATP_FBAsolutionM.x(BiomassMISO));
ISO28 = findRxnIDs(modelM_Ratio,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', ATP_FBAsolutionM.x(ISO28));
ISO40 = findRxnIDs(modelM_Ratio,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', ATP_FBAsolutionM.x(ISO40));
ISO60 = findRxnIDs(modelM_Ratio,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', ATP_FBAsolutionM.x(ISO60));
ISOM = findRxnIDs(modelM_Ratio,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', ATP_FBAsolutionM.x(ISOM));



