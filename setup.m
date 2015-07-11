%% Set paths
addpath('external/MinMaxSelection');
addpath('util');

imgset  = [class, '_set_' num2str(setid)];
imgDir  = [baseDir, '/data/' imgset '/'];
resDir  = [baseDir '/results/' imgset '/flowweb/'];
initDir = [baseDir '/results/' imgset '/dsp/'];

if ~exist(resDir, 'dir')
    mkdir(resDir);
end
if gpu > 0
    g = gpuDevice(gpu);
end

%% Set parameters
params = [];
% Relaxed threshold (normalized w.r.t. image size) for cycle completion
params.cycleThresh = 0.05;
% FIXME: might need to tune the following
% Only top ranked flows from each pairwise flow field is allowed to enter
% the overall priority sorting. This is mainly for speed up and memory 
% efficiency
params.proposalRate = 0.2;
% Determines how many flows are updated during inter-image phase (referred 
% to as \beta\% in the paper)
params.updateRate = 0.2;
% Number of maximum iterations allowed
params.maxIter = 10;
% Terminate the algorithm if cycle improvement is below the threshold
params.convergeThresh = 0.03;
% How much spatial proximity affects in the filtering/intra-image phase
params.filterSpatialSigma = 0.05;
% How much cycle consistency affects in the filtering/intra-image phase
params.filterCycleSigma = 0.05;
% Determines how many flows are filtered during intra-image phase. Only
% flows of relatively low cycle-consistency are filtered
params.filterRate = 0.5;
% Threshold for non-zero spatial weights in the filtering/intra-image phase
params.spatialWeightThresh = 0.0001;
% Regularization towards the initial solution
params.lambda = 0.0;
params.resDir = resDir;
params.saveEveryIter = true;