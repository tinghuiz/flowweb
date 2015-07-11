% metric = 'keypoint';
% class = 'car';
% setid = 2;
% method = 'dsp';
baseDir = 'pascal/';
imgset  = [class, '_set_' num2str(setid)];
annoDir = [baseDir, '/data/' imgset '/'];
annoFiles = dir([annoDir, '*.mat']);
annos = cell(1, length(annoFiles));
addpath('util');
for i = 1 : length(annoFiles)
    load([annoDir, annoFiles(i).name]);
    if strcmp(metric, 'keypoint')
        annos{i}.keypts = keypts;
        annos{i}.keyptStatus = keypts_status;
    end
    if strcmp(metric, 'part')
        annos{i}.parts = parts;
    end
end
flowDir  = [baseDir '/results/' imgset '/' method '/'];
flowFile = [flowDir, 'pairflows.mat'];
load(flowFile);

params.metric = metric;
params.alpha = 0.05;
params.numImage = size(pairvx,1);
params.height = size(pairvx{1,2},1);
params.width = size(pairvx{1,2},2);
res = evaluate(pairvx, pairvy, annos, params);
fprintf('class = %s, setid = %d, method = %s, metric = %s, res = %f\n', ...
    class, setid, method, metric, res);