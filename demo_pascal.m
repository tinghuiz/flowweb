class = 'car';
setid = 2;
baseDir = 'pascal/';

setup;

ims = loadImageSet(imgDir);
[pairvx, pairvy] = initflows(ims, initDir);
[pairvx, pairvy] = alignFlowWeb(pairvx, pairvy, params);
save([resDir, 'pairflows.mat'], 'pairvx', 'pairvy');
