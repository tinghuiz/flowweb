%% Set paths
addpath('external/MinMaxSelection');
addpath('external/pwmetric');
addpath('util');

imgset  = [class, '_set_' num2str(setid)];
imgDir  = [baseDir, '/data/' imgset '/'];
resDir  = [baseDir '/results/' imgset '/flowweb/'];
initDir = [baseDir '/results/' imgset '/' initMethod '/'];

if ~exist(resDir, 'dir')
    mkdir(resDir);
end

%% Set parameters
params = [];
% Relaxed threshold (normalized w.r.t. image size) for cycle completion
params.cycleThresh = 0.05;
% Only top ranked flows from each pairwise flow field is allowed to enter
% the overall priority sorting. For speed-up purpose.
params.proposalRate = 0.4;
% Determines how many flows are updated during inter-image phase
params.updateRate = 0.2;
% Number of maximum iterations allowed
params.maxIter = 10;
% Terminate the algorithm if cycle improvement is below the threshold
params.convergeThresh = 0.03;
% Determines the spatial weight in the intra-image phase
params.filterSpatialSigma = 0.05;
% Determines the consistency weight in the intra-image phase
params.filterCycleSigma = 0.05;
% Determines how many flows are filtered during intra-image phase. Only
% flows of relatively low cycle-consistency are filtered
params.filterRate = 0.25;
% Regularization towards the initial solution
params.lambda = 0.01;
% Whether to save the result of each iteration
params.saveEveryIter = false;
params.resDir = resDir;

if gpu > 0
    params.isGPU = true;
    g = gpuDevice(gpu);
else
    params.isGPU = false;
end