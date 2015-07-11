%% Convert the old mat results to cell
metric = 'part';
class = 'car';
setid = 2;
method = 'fweb';
baseDir = 'pascal/';
imgset  = [class, '_set_' num2str(setid)];
imgDir  = [baseDir, '/data/' imgset '/'];
flowDir  = [baseDir '/results/' imgset '/' method '/'];
flowFile = [flowDir, 'pairflows.mat'];
load(flowFile);
im = imread([imgDir, '001.jpg']);
[pairvx, pairvy] = flowMatToCell(pairvx, pairvy, size(im,1));
save(flowFile, 'pairvx', 'pairvy');