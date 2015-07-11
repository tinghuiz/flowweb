function [pairvx, pairvy] = cyclefilter(pairvx, pairvy, cycleSet, ...
    nnInds, nnW, initvx, initvy, height, width, params)

N = sqrt(size(pairvx, 2));
cycleMap = int16(sum(cycleSet, 3));

for src = 1 : N
    if src > 1
        fprintf('%.3d/%.3d', src, N);
    end
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        s2t = src + (tgt - 1) * N;
        nnCycle = reshape(cycleMap(nnInds(:), s2t), size(nnInds));
        cycleDiff = bsxfun(@minus, nnCycle, cycleMap(:, s2t));
        cycleWeights = exp(single(cycleDiff)/...
            params.filterCycleSigma/(N-2));
        nnvx = reshape(pairvx(nnInds(:), s2t), size(nnInds));
        nnvy = reshape(pairvy(nnInds(:), s2t), size(nnInds));
        nnReg = fdist(nnvx, nnvy, repmat(initvx(:, s2t), 1, ...
            size(nnInds,2)), repmat(initvy(:, s2t), 1, size(nnInds,2)), ...
            height, width);
        directReg = fdist(pairvx(:, s2t), pairvy(:, s2t),...
            initvx(:, s2t), initvy(:, s2t), height, width);
        regWeights = params.lambda * (nnReg - repmat(directReg, ...
            1, size(nnInds,2)));
        cycleWeights = cycleWeights - regWeights;
        cycleWeights(cycleWeights < 1) = 0;
        weightsAll = nnW .* cycleWeights;
        weightsAll = bsxfun(@times, weightsAll, 1./sum(weightsAll, 2));
        nnvx = reshape(pairvx(nnInds(:), s2t), size(nnInds));
        pairvx(:, s2t) = sum(nnvx .* weightsAll, 2);
        nnvy = reshape(pairvy(nnInds(:), s2t), size(nnInds));
        pairvy(:, s2t) = sum(nnvy .* weightsAll, 2);
    end
    if src > 1
        % FIXME: get number of backspace needed automatically
        fprintf(repmat('\b',1,7));
    end
end
pairvx = round(pairvx);
pairvy = round(pairvy);