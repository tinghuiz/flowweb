% class = 'car';
% setid = 100;
baseDir = 'pascal/';
gpu = 1;

setup;
ims = loadImageSet(imgDir);
% Initialize FlowWeb
[pairvx, pairvy] = initflows(ims, initDir);
% Joint alignment
[pairvx, pairvy] = alignFlowWeb(pairvx, pairvy, params);
% Save aligned FlowWeb
% Use int16 for space saving
pairvx = int16(round(pairvx));
pairvy = int16(round(pairvy));
[pairvx, pairvy] = flowMatToCell(pairvx, pairvy, size(ims{1},1));
save([resDir, 'pairflows.mat'], 'pairvx', 'pairvy');