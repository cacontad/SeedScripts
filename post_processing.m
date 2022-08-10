initCobraToolbox
% cd to working directory
% perform sampling analysis 
warmupn = 6070; % min required according to our setting
nFiles = 170; % it can be modified according your requirements
pointsPerFile = 5000;
stepsPerPoint = 500;
fileBaseNo = 0;
maxTime = 3600000;

% each conditional model should be used to run the sampling
load('Sampling_C08_XX.mat', 'model_X')
path ='PATH to save your sampling files';
fileName = 'Sampling_X';
performSampling(model_X, warmupn, fileName, nFiles, pointsPerFile, stepsPerPoint, fileBaseNo, maxTime, path);

% or if sampling files are already available, load files
% number of Files generated during sampling
nFiles = 170;
pointsPerFile = 5000;
path = '/PATH to your sampling files';
% root name of your sampling files
fileName = '/Sampling_XX';
samples = loadSamples([path fileName], nFiles, pointsPerFile);

% Load your condition specific model
load('Sampling_C08_XX.mat', 'model_X')

% setting folder for post-processing
fileBaseNo = 0;
fileName= 'fluxWoCycle_C08';
path ='/PATH to save post-processing';
fileName = [path filesep fileName];

% update according your number of sampling files
a = zeros(4,170);
% update model according to the sampling data
model_A = model_X;

for i = 1:170

testing_nocycles = samples(:,(i-1)*5000+1:(i-1)*5000+5000);

% feasTol can be ajusted
feasTol = 10 .^ (-7:-4);
fixRxns = model_A.rxns(model_A.lb>0);

K = size(testing_nocycles, 2);
rxnEx = sum(model_A.S ~= 0, 1) <= 1;
checkRxnIDs = [find(rxnEx(:)); findRxnIDs(model_A, fixRxns(:))];

SUCCESSs = false(K, numel(feasTol));
success = false(K, 1);
id = 1:K;
fluxWoCycle = zeros(size(testing_nocycles));


    for j = 1:4
        if sum(a(:,i),1) >= 5000
        break
        end
        id = id(~success);
        [fluxWoCycleJ, success] = removeCyclesCobra(model_A, testing_nocycles(:, id), 0, fixRxns, feasTol(j));
        SUCCESSs(id(success), j) = true;
        dev = testing_nocycles(checkRxnIDs, id(success)) - fluxWoCycleJ(checkRxnIDs, success);
        infeas1 = checkSolFeas(model_A, testing_nocycles(:, id(success)));
        infeas2 = checkSolFeas(model_A, fluxWoCycleJ(:, success));
        fluxWoCycle(:, id(success)) = fluxWoCycleJ(:, success);
        a(j,i) = sum(success);
        fprintf('feasTol: %.4e\n', feasTol(j))
        fprintf('Max deviation in exchange fluxes and fixed internal fluxes: %.6e\n', max(max(dev)))
        fprintf('Max infeasibility in the original flux distributions: %.6e\n', max(infeas1))
        fprintf('Max infeasibility in the new flux distributions: %.6e\n\n', max(infeas2))
        
    end
points = fluxWoCycle;  
file = [fileName '_' num2str(fileBaseNo+i) '.mat'];
save (file,'points');
end

clearvars -EXCEPT a
save Success_Ratio_1_1_MP_hetero_C08_M.mat a;
