% ATP is enforced
% No Maternal contribution
clear 
pathModel = 'YOUR_PATH';
filename = 'C08_seed.mat';
load([pathModel, filename])
model = C08_seed;
clear C08_seed

%% ISO MAX
%% 28 DAF
model28 = model;

% ATPM can be adjusted according to ATP:NADPH ratios
ATPMaintenance = find(ismember(model.rxns,'RXN-11109_c'));
model28.lb(ATPMaintenance) = 3.1662;

Biomass = find(ismember(model28.rxns,'BIOMASS-RXN_C08_28DAF'));
GrowthPerHour = 0.306275525/24;
model28.lb(Biomass) = GrowthPerHour;
model28.ub(Biomass) = GrowthPerHour;


%% Nutrients according 28DAF: literature data
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

%% One Active Biomass rxn

% rxnID = findRxnIDs(model4,'BIOMASS-RXN_C08_28DAF');
% model4.lb(rxnID,1) = 0;
% model4.ub(rxnID,1) = 0;

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

%% Iso out of the model
model28 = addReaction(model28,{'sink_Isoflavones_28DAF[c]','Isoflavones_28DAF[c]'},{'Isoflavones_28DAF[c]'},-1,false);
model28 = changeObjective(model28,'sink_Isoflavones_28DAF[c]',1);

%% OF
model28.osenseStr = char('max');
C08_ISO_FBAsolution28 = optimizeCbModel(model28);


fprintf('\nC08-max ISO_28 flux through the network : %.10f\n', C08_ISO_FBAsolution28.f);
fprintf('28 DAF Biomass: %.4f\n', C08_ISO_FBAsolution28.x(Biomass));
fprintf('ATPMaintenance: %.4f\n', C08_ISO_FBAsolution28.x(ATPMaintenance));


ISO28 = findRxnIDs(model28,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', C08_ISO_FBAsolution28.x(ISO28));
ISO40 = findRxnIDs(model28,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', C08_ISO_FBAsolution28.x(ISO40));
ISO60 = findRxnIDs(model28,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', C08_ISO_FBAsolution28.x(ISO60));
ISOM = findRxnIDs(model28,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', C08_ISO_FBAsolution28.x(ISOM));
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model28, C08_ISO_FBAsolution28.w, 1, 1)

model28.lb(ATPMaintenance) = 0;
ISO_FBAsolution28 = optimizeCbModel(model28);
fprintf('\nC08-max ISO_28 flux without ATPM: %.10f\n', ISO_FBAsolution28.f);
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model28, ISO_FBAsolution28.w, 1, 1)
fprintf('\n');

%% 40DAF
model40 = model;

% ATPM can be adjusted according to ATP:NADPH ratios
ATPMaintenance = find(ismember(model.rxns,'RXN-11109_c'));
model40.lb(ATPMaintenance) = 2.9322;

Biomass = find(ismember(model40.rxns,'BIOMASS-RXN_C08_40DAF'));
GrowthPerHour = 0.471988925/24;
model40.lb(Biomass) = GrowthPerHour;
model40.ub(Biomass) = GrowthPerHour;

%% Nutrients according 40DAF:
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

%% One Active Biomass rxn

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

%% Iso out of the model
model40 = addReaction(model40,{'sink_Isoflavones_40DAF[c]','Isoflavones_40DAF[c]'},{'Isoflavones_40DAF[c]'},-1,false);
model40 = changeObjective(model40,'sink_Isoflavones_40DAF[c]',1);

%% OF
model40.osenseStr = char('max');
C08_ISO_FBAsolution40 = optimizeCbModel(model40);


fprintf('\nC08-max ISO_40 flux through the network : %.10f\n', C08_ISO_FBAsolution40.f);
fprintf('40 DAF Biomass: %.4f\n', C08_ISO_FBAsolution40.x(Biomass));
fprintf('ATPMaintenance: %.4f\n', C08_ISO_FBAsolution40.x(ATPMaintenance));

ISO28 = findRxnIDs(model40,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', C08_ISO_FBAsolution40.x(ISO28));
ISO40 = findRxnIDs(model40,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', C08_ISO_FBAsolution40.x(ISO40));
ISO60 = findRxnIDs(model40,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', C08_ISO_FBAsolution40.x(ISO60));
ISOM = findRxnIDs(model40,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', C08_ISO_FBAsolution40.x(ISOM));
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model40, C08_ISO_FBAsolution40.w, 1, 1)

model40.lb(ATPMaintenance) = 0;
ISO_FBAsolution40 = optimizeCbModel(model40);
fprintf('\nC08-max ISO_40 flux without ATPM: %.10f\n', ISO_FBAsolution40.f);
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model40, ISO_FBAsolution40.w, 1, 1)
fprintf('\n');
%% 60DAF

model60 = model;

% ATPM can be adjusted according to ATP:NADPH ratios
ATPMaintenance = find(ismember(model.rxns,'RXN-11109_c'));
model60.lb(ATPMaintenance) = 0.0506;

Biomass = find(ismember(model60.rxns,'BIOMASS-RXN_C08_60DAF'));
GrowthPerHour = 0.258666711/24;
model60.lb(Biomass) = GrowthPerHour;
model60.ub(Biomass) = GrowthPerHour;

%% Nutrients according 60DAF:
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
model60.lb(rxnID,1) = 0;

%% One Active Biomass rxn

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

%% Iso out of the model
model60 = addReaction(model60,{'sink_Isoflavones_60DAF[c]','Isoflavones_60DAF[c]'},{'Isoflavones_60DAF[c]'},-1,false);
model60 = changeObjective(model60,'sink_Isoflavones_60DAF[c]',1);

%% OF
model60.osenseStr = char('max');
C08_ISO_FBAsolution60 = optimizeCbModel(model60);

fprintf('\nC08-max ISO_60 flux through the network : %.10f\n', C08_ISO_FBAsolution60.f);
fprintf('60 DAF Biomass: %.4f\n', C08_ISO_FBAsolution60.x(Biomass));
fprintf('ATPMaintenance: %.4f\n', C08_ISO_FBAsolution60.x(ATPMaintenance));

ISO28 = findRxnIDs(model60,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', C08_ISO_FBAsolution60.x(ISO28));
ISO40 = findRxnIDs(model60,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', C08_ISO_FBAsolution60.x(ISO40));
ISO60 = findRxnIDs(model60,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', C08_ISO_FBAsolution60.x(ISO60));
ISOM = findRxnIDs(model60,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', C08_ISO_FBAsolution60.x(ISOM));
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model60, C08_ISO_FBAsolution60.w, 1, 1)

model60.lb(ATPMaintenance) = 0;
ISO_FBAsolution60 = optimizeCbModel(model60);
fprintf('\nC08-max ISO_60 flux without ATPM: %.10f\n', ISO_FBAsolution60.f);
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model60, ISO_FBAsolution60.w, 1, 1)
fprintf('\n');
%% MATURE

modelM = model;

% ATPM can be adjusted according to ATP:NADPH ratios
ATPMaintenance = find(ismember(model.rxns,'RXN-11109_c'));
modelM.lb(ATPMaintenance) = 0.0819;

Biomass = find(ismember(modelM.rxns,'BIOMASS-RXN_C08_MATURE'));
GrowthPerHour = 0.223810871/24;
modelM.lb(Biomass) = GrowthPerHour;
modelM.ub(Biomass) = GrowthPerHour;

%% Nutrients according MATURE:
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
modelM.lb(rxnID,1) = 0;

%% One Active Biomass rxn

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

%% Iso out of the model
modelM = addReaction(modelM,{'sink_Isoflavones_MATURE[c]','Isoflavones_MATURE[c]'},{'Isoflavones_MATURE[c]'},-1,false);
modelM = changeObjective(modelM,'sink_Isoflavones_MATURE[c]',1);

%% OF
modelM.osenseStr = char('max');
C08_ISO_FBAsolutionM = optimizeCbModel(modelM);



fprintf('\nC08-max ISO_M flux through the network : %.10f\n', C08_ISO_FBAsolutionM.f);
fprintf('MATURE Biomass: %.4f\n', C08_ISO_FBAsolutionM.x(Biomass));
fprintf('ATPMaintenance: %.4f\n', C08_ISO_FBAsolutionM.x(ATPMaintenance));

ISO28 = findRxnIDs(modelM,'ISOFLAVONES-RXN_28DAF');
fprintf('ISOFLAVONES-RXN_28DAF flux: %.10f\n', C08_ISO_FBAsolutionM.x(ISO28));
ISO40 = findRxnIDs(modelM,'ISOFLAVONES-RXN_40DAF');
fprintf('ISOFLAVONES-RXN_40DAF flux: %.10f\n', C08_ISO_FBAsolutionM.x(ISO40));
ISO60 = findRxnIDs(modelM,'ISOFLAVONES-RXN_60DAF');
fprintf('ISOFLAVONES-RXN_60DAF flux: %.10f\n', C08_ISO_FBAsolutionM.x(ISO60));
ISOM = findRxnIDs(modelM,'ISOFLAVONES-RXN_MATURE');
fprintf('ISOFLAVONES-RXN_MATURE flux: %.10f\n', C08_ISO_FBAsolutionM.x(ISOM));
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(modelM, C08_ISO_FBAsolutionM.w, 1, 1)



modelM.lb(ATPMaintenance) = 0;
ISO_FBAsolutionM = optimizeCbModel(modelM);
fprintf('\nC08-max ISO_M flux without ATPM: %.10f\n', ISO_FBAsolutionM.f);
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(modelM, ISO_FBAsolutionM.w, 1, 1)

