class = 'car';
setid = 1;
baseDir = 'pascal_toy/';
initMethod = 'dsp';
% The intra-image phase is GPU accelerated
gpu = 1;
% Setup relavent paths and parameters
alignSetup;
% Load the input image set
ims = loadImageSet(imgDir);
% Initialize FlowWeb
[pairvx, pairvy] = initflows(ims, initDir);
% Joint alignment
[pairvx, pairvy] = jointAlign(pairvx, pairvy, params);
% Save aligned FlowWeb
[pairvx, pairvy] = flowMatToCell(pairvx, pairvy, size(ims{1},1));
save([resDir, 'pairflows.mat'], 'pairvx', 'pairvy');

% Evaluation
metric = 'keypoint'; % can also be 'part'

evalParams = [];
evalParams.metric = metric;
evalParams.alpha = 0.05;
evalParams.numImage = length(ims);
evalParams.height = size(ims{1},1);
evalParams.width = size(ims{1},2);
evalParams.resRoot = [baseDir '/results/' imgset '/'];
evalParams.annoDir = imgDir;

evalParams.method = 'dsp';
resInit = evaluate(evalParams);
evalParams.method = 'flowweb';
resOurs = evaluate(evalParams);
fprintf('class = %s, setid = %d,  metric = %s\n', class, setid, metric);
fprintf('DSP (initialization) = %.3f, Ours = %.3f\n', resInit, resOurs)

