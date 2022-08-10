clear

pathModel = 'YOUR_PATH';
filename = 'W05_seed.mat';
load([pathModel, filename])
model = W05_seed;
clear W05_seed

model.osenseStr = char('max');
model = changeObjective(model,'BIOMASS-RXN_W05_28DAF',1);

%% no biomass without nutrients
model2 = model;
[selExc,selUpt] = findExcRxns(model2);

Ex_ind = find(selExc);
for i = 1:length(Ex_ind)
    RxnInd = Ex_ind(i);
    model2.lb(RxnInd) = 0;
end

% 28DAF biomass rxn
model28 = model2;
model28 = changeObjective(model28,'BIOMASS-RXN_W05_28DAF',1);
s = optimizeCbModel(model28);
fprintf('\n28 DAF WT growth rate no nutrients: %.4f\n', s.f)
model28 = changeObjective(model28,'BIOMASS-RXN_W05_28DAF_ISOFLAVONES',1);
s = optimizeCbModel(model28);
fprintf('\n28 DAF_ISO WT growth rate no nutrients: %.4f\n', s.f)

% 40DAF biomass rxn
model40 = model2;
model40 = changeObjective(model40,'BIOMASS-RXN_W05_40DAF',1);
s = optimizeCbModel(model40);
fprintf('\n40 DAF WT growth rate no nutrients: %.4f\n', s.f)
model40 = changeObjective(model40,'BIOMASS-RXN_W05_40DAF_ISOFLAVONES',1);
s = optimizeCbModel(model40);
fprintf('\n40 DAF_ISO WT growth rate no nutrients: %.4f\n', s.f)

% 60DAF biomass rxn
model60 = model2;
model60 = changeObjective(model60,'BIOMASS-RXN_W05_60DAF',1);
s = optimizeCbModel(model60);
fprintf('\n60 DAF WT growth rate no nutrients: %.4f\n', s.f)
model60 = changeObjective(model60,'BIOMASS-RXN_W05_60DAF_ISOFLAVONES',1);
s = optimizeCbModel(model60);
fprintf('\n60 DAF_ISO WT growth rate no nutrients: %.4f\n', s.f)

% MATURE biomass rxn
modelM = model2;
modelM = changeObjective(modelM,'BIOMASS-RXN_W05_MATURE',1);
s = optimizeCbModel(modelM);
fprintf('\nMATURE WT growth rate no nutrients: %.4f\n', s.f)
modelM = changeObjective(modelM,'BIOMASS-RXN_W05_MATURE_ISOFLAVONES',1);
s = optimizeCbModel(modelM);
fprintf('\nMATURE_ISO WT growth rate no nutrients: %.4f\n', s.f)

%% biomass production with unlimited nutrients

% 28DAF biomass rxn
model28 = model;
model28 = changeObjective(model28,'BIOMASS-RXN_W05_28DAF',1);
s = optimizeCbModel(model28,'max','one');
fprintf('\n 28DAF Biomass rxn WT growth rate unlimited nutrients: %.4f\n', s.f)
model28 = changeObjective(model28,'BIOMASS-RXN_W05_28DAF_ISOFLAVONES',1);
s = optimizeCbModel(model28);
fprintf('\n28 DAF_ISO WT growth rate unlimited nutrients: %.4f\n', s.f)

% 40DAF biomass rxn
model40 = model;
model40 = changeObjective(model40,'BIOMASS-RXN_W05_40DAF',1);
s = optimizeCbModel(model40,'max','one');
fprintf('\n 40DAF Biomass rxn WT growth rate unlimited nutrients: %.4f\n', s.f)
model40 = changeObjective(model40,'BIOMASS-RXN_W05_40DAF_ISOFLAVONES',1);
s = optimizeCbModel(model40);
fprintf('\n40 DAF_ISO WT growth rate unlimited nutrients: %.4f\n', s.f)

% 60DAF biomass rxn
model60 = model;
model60 = changeObjective(model60,'BIOMASS-RXN_W05_60DAF',1);
s = optimizeCbModel(model60,'max','one');
fprintf('\n 60DAF Biomass rxn WT growth rate unlimited nutrients: %.4f\n', s.f)
model60 = changeObjective(model60,'BIOMASS-RXN_W05_60DAF_ISOFLAVONES',1);
s = optimizeCbModel(model60);
fprintf('\n60 DAF_ISO WT growth rate unlimited nutrients: %.4f\n', s.f)

% Mature biomass rxn
modelM = model;
modelM = changeObjective(modelM,'BIOMASS-RXN_W05_MATURE',1);
s = optimizeCbModel(modelM,'max','one');
fprintf('\n Mature Biomass rxn WT growth rate unlimited nutrients: %.4f\n', s.f)
modelM = changeObjective(modelM,'BIOMASS-RXN_W05_MATURE_ISOFLAVONES',1);
s = optimizeCbModel(modelM);
fprintf('\nMATURE_ISO WT growth rate unlimited nutrients: %.4f\n', s.f)

%% limited nutrients
model4 = model;

% 28 DAF conditions: from Initial_setup excel file
% literature data: [umol g DW-1 h-1]
rxnID = findRxnIDs(model,'Exchange_Sucrose');
model4.lb(rxnID,1) = -41;
rxnID = findRxnIDs(model,'Exchange_GLC');
model4.lb(rxnID,1) = -5;
rxnID = findRxnIDs(model,'Exchange_GLN');
model4.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model,'Exchange_ASN');
model4.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model,'Exchange_O2');
model4.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model,'Exchange_Photon');
model4.lb(rxnID,1) = -10;

s3 = optimizeCbModel(model4);
fprintf('\nWT growth rate 28DAF limited nutrients: %.4f\n', s3.f)
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model4, s3.w, 1, 1)

% no glutamine or asparagine
model4_1 = changeRxnBounds(model4, {'Exchange_GLN','Exchange_ASN'}, 0, 'l');
s1 = optimizeCbModel(model4_1);
fprintf('\nWT growth rate 28DAF without glutamine and asparagine: %.4f\n', s1.f)
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model4_1, s1.w, 1, 1)

% no sucrose or glucose
model4_2 = changeRxnBounds(model4, {'Exchange_Sucrose','Exchange_GLC'}, 0, 'l');
s2 = optimizeCbModel(model4_2);
fprintf('\nWT growth rate 28DAF without sucrose and glucose: %.4f\n', s2.f)
fprintf('Nonzero shadow price (limiting substrates):\n');
printFluxVector(model4_2, s2.w, 1, 1)

%% Isoflavonoid production
model5 = model;

% 28 DAF conditions
rxnID = findRxnIDs(model,'Exchange_Sucrose');
model5.lb(rxnID,1) = -41;
rxnID = findRxnIDs(model,'Exchange_GLC');
model5.lb(rxnID,1) = -5;
rxnID = findRxnIDs(model,'Exchange_GLN');
model5.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model,'Exchange_ASN');
model5.lb(rxnID,1) = -10;
rxnID = findRxnIDs(model,'Exchange_O2');
model5.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model,'Exchange_Photon');
model5.lb(rxnID,1) = -10;

% CPD-3421[c]	genistin
model5 = addReaction(model5,{'sink_CPD-3421[c]','sink_CPD-3421[c]'},{'CPD-3421[c]'},-1,false);
model5 = changeObjective(model5,'sink_CPD-3421[c]',1);
s = optimizeCbModel(model5,'max','one');
fprintf('\nGenistin 28DAF limited nutrients: %.4f\n', s.f)

% CPD-3424[c]	daidzin
model5 = addReaction(model5,{'sink_CPD-3424[c]','sink_CPD-3424[c]'},{'CPD-3424[c]'},-1,false);
model5 = changeObjective(model5,'sink_CPD-3424[c]',1);
s = optimizeCbModel(model5,'max','one');
fprintf('\nDaidzin 28DAF limited nutrients: %.4f\n', s.f)

% GLYCITIN[c]
model5 = addReaction(model5,{'sink_GLYCITIN[c]','sink_GLYCITIN[c]'},{'GLYCITIN[c]'},-1,false);
model5 = changeObjective(model5,'sink_GLYCITIN[c]',1);
s = optimizeCbModel(model5,'max','one');
fprintf('\nGlycitin 28DAF limited nutrients: %.4f\n', s.f)

%% Day and night conditions
model6 = model;

% Heterotrophic condition (NIGHT)
% From Soybean paper:
% "The heterotrophic condition was modeled by setting photon influx to zero 
% and allowing Glc as the carbon and energy input (GLC_tx), with all other 
% constraints the same as in the phototrophic condition."
rxnID = findRxnIDs(model,'Exchange_Sucrose');
model6.lb(rxnID,1) = 0;
rxnID = findRxnIDs(model,'Exchange_GLC');
model6.lb(rxnID,1) = -1000;
rxnID = findRxnIDs(model,'Exchange_GLN');
model6.lb(rxnID,1) = 0;
rxnID = findRxnIDs(model,'Exchange_ASN');
model6.lb(rxnID,1) = 0;
rxnID = findRxnIDs(model,'Exchange_Photon');
model6.lb(rxnID,1) = 0;

s = optimizeCbModel(model6,'max','one');
fprintf('\nNight: %.4f\n', s.f)

printFluxVector(model6, s.x, 1, 1)

rxnID = findRxnIDs(model6,'RXN490-3650_p');
fprintf('\nFlux photosystem I: %.4f\n', s.x(rxnID,1))
rxnID = findRxnIDs(model6,'PSII-RXN_p');
fprintf('\nFlux photosystem II: %.4f\n', s.x(rxnID,1))
rxnID = findRxnIDs(model6,'Photon_tx');
fprintf('\nFlux Photon Uptake: %.4f\n', s.x(rxnID,1))

% Phototrophic condition (DAY)
% From Soybean paper:
% "In phototrophic conditions, a photon influx was allowed from 
% the photon input reaction (Photon_tx) with active photorespiration 
% and with all organic sources of energy set to zero. 
% Mineral nutrients were allowed to exchange with the environment."
rxnID = findRxnIDs(model,'Exchange_GLC');
model6.lb(rxnID,1) = 0;
rxnID = findRxnIDs(model,'Exchange_Photon');
model6.lb(rxnID,1) = -1000;
rxnID = findRxnIDs(model,'Exchange_CO2');
model6.lb(rxnID,1) = -1000;


s = optimizeCbModel(model6,'max','one');
fprintf('\nDay: %.4f\n', s.f)


printFluxVector(model6, s.x, 1, 1)

rxnID = findRxnIDs(model6,'RXN490-3650_p');
fprintf('\nFlux photosystem I: %.4f\n', s.x(rxnID,1))
rxnID = findRxnIDs(model6,'PSII-RXN_p');
fprintf('\nFlux photosystem II: %.4f\n', s.x(rxnID,1))
rxnID = findRxnIDs(model6,'Photon_tx');
fprintf('\nFlux Photon Uptake: %.4f\n', s.x(rxnID,1))

%% Mature metabolites
model10 = model;

% Mature conditions: [umol g DW-1 h-1]
rxnID = findRxnIDs(model10,'Exchange_Sucrose');
model10.lb(rxnID,1) = -15;
rxnID = findRxnIDs(model10,'Exchange_GLC');
model10.lb(rxnID,1) = -2;
rxnID = findRxnIDs(model10,'Exchange_GLN');
model10.lb(rxnID,1) = -17;
rxnID = findRxnIDs(model10,'Exchange_ASN');
model10.lb(rxnID,1) = -17;
rxnID = findRxnIDs(model10,'Exchange_O2');
model10.lb(rxnID,1) = -1;
rxnID = findRxnIDs(model10,'Exchange_Photon');
model10.lb(rxnID,1) = -1;

fid = fopen('matureMetabolites.txt');
tline = fgetl(fid);
i = 0;
while ischar(tline)
    i = i+1;
    metaboliteMature(i,1) = cellstr(tline);
    tline = fgetl(fid);
end
fclose(fid);

for j = 1:numel(metaboliteMature)
    model10 = addReaction(model10,['Sink_' metaboliteMature{j,1}], 'reactionFormula', [metaboliteMature{j,1}, ' <=>']);
end

rxnSink = strcat('Sink_', metaboliteMature);
gr = zeros(numel(rxnSink), 2);

for j = 1:numel(metaboliteMature)
    model10 = changeRxnBounds(model10, rxnSink, 0, 'b');
    model10 = changeRxnBounds(model10, rxnSink{j}, 1000, 'u');
    model10 = changeObjective(model10,rxnSink{j},1);
    s = optimizeCbModel(model10);
    gr(j, 2) = s.f;
    fprintf('Production of %s: %.2f\n', metaboliteMature{j}, gr(j,2));   
end

