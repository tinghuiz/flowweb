class = 'car';
setid = 100;
baseDir = 'pascal/';

setup;

ims = loadImageSet(imgDir);
[pairvx, pairvy] = initflows(ims, initDir);
[pairvx, pairvy] = alignFlowWeb(pairvx, pairvy, params);
save([resDir, 'pairflows.mat'], 'pairvx', 'pairvy');
