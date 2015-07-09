function [pairvx, pairvy] = initflows(ims, initDir)

N = length(ims);

if exist([initDir, 'pairflows.mat'], 'file')
    load([initDir, 'pairflows.mat']);
else
    % By default use DSP to initialize pairwise flows
    addpath('external/dsp-code/');
    vlfeat_root = 'external/dsp-code/vlfeat-0.9.17/';
    try
        vl_root();
    catch
        run([vlfeat_root, 'toolbox/vl_setup']);
    end
    if ~exist(initDir, 'dir')
        mkdir(initDir);
    end
    pairvx = cell(N,N);
    pairvy = cell(N,N);
    feats = cell(1, N);
    for i = 1 : N
        feats{i} = ExtractSIFT_WithPadding(ims{i}, [], 4);
    end
    for i = 1 : N
        fprintf('Init progress: %d/%d\n', i, N);
        parfor j = 1 : N
            if i == j
                continue;
            end
            [pairvx{i,j}, pairvy{i,j}] = DSPMatch(feats{i}, feats{j});
            pairvx{i,j} = single(pairvx{i,j});
            pairvy{i,j} = single(pairvy{i,j});
        end
    end
    save([initDir, 'pairflows.mat'], 'pairvx', 'pairvy');
end