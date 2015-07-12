function [pairvx, pairvy] = cyclefilter(pairvx, pairvy, cycleSet, ...
    spatialWeights, params, is_gpu)

N = params.numImage;
% Per-flow consistency score
cycleMap = int16(sum(cycleSet, 3));
% Number of flows to be filtered 
numfilter = round(params.filterRate * params.width * params.height);

if is_gpu
    cycleMap = gpuArray(cycleMap);
    spatialWeights = gpuArray(spatialWeights);
    pairvx = gpuArray(pairvx);
    pairvy = gpuArray(pairvy);
else
    pairvx = gather(pairvx);
    pairvy = gather(pairvy);
    spatialWeights = gather(spatialWeights);
end

for src = 1 : N
    if src > 1
        fprintf('%.3d/%.3d', src, N);
    end
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        s2t = src + (tgt - 1) * N;
        [~, fidx] = mink(gather(cycleMap(:, s2t)), numfilter);
        cycleDiff = bsxfun(@minus, cycleMap(:, s2t)', ...
           cycleMap(fidx, s2t));
        cycleWeights = exp(single(cycleDiff)/params.filterCycleSigma/(N-2));
        cycleWeights(cycleWeights < 1) = 0;
        weightsAll = spatialWeights(fidx, :) .* cycleWeights;
        weightsAll = bsxfun(@times, weightsAll, 1./sum(weightsAll, 2));
        pairvx(fidx, s2t) = sum(bsxfun(@times, weightsAll, pairvx(:, s2t)'),2);
        pairvy(fidx, s2t) = sum(bsxfun(@times, weightsAll, pairvy(:, s2t)'),2);
    end
    if src > 1
        % FIXME: get number of backspace needed automatically
        fprintf(repmat('\b',1,7));
    end
end
pairvx = round(pairvx);
pairvy = round(pairvy);
pairvx = gather(pairvx);
pairvy = gather(pairvy);

% function [pairvx, pairvy] = cyclefilter(pairvx, pairvy, cycleSet, ...
%     nbInds, nnW, initvx, initvy, height, width, params)
% 
% N = sqrt(size(pairvx, 2));
% cycleMap = int16(sum(cycleSet, 3));
% numfilter = round(params.filterRate * params.width * params.height);
% 
% for src = 1 : N
%     if src > 1
%         fprintf('%.3d/%.3d', src, N);
%     end
%     for tgt = 1 : N
%         if src == tgt
%             continue;
%         end
%         s2t = src + (tgt - 1) * N;
%         [~, thisInds] = mink(cycleMap(:, s2t), numfilter);
%         nb = nbInds(thisInds, :);
%         nbCycle = reshape(cycleMap(nb(:), s2t), size(nb));
%         cycleDiff = bsxfun(@minus, nbCycle, cycleMap(thisInds, s2t));
%         cycleWeights = exp(single(cycleDiff)/...
%             params.filterCycleSigma/(N-2));
%         nbvx = reshape(pairvx(nb(:), s2t), size(nb));
%         nbvy = reshape(pairvy(nb(:), s2t), size(nb));
%         nbReg = fdist(nbvx, nbvy, repmat(initvx(thisInds, s2t), 1, ...
%             size(nb,2)), repmat(initvy(thisInds, s2t), 1, size(nb,2)), ...
%             height, width);
%         directReg = fdist(pairvx(thisInds, s2t), pairvy(thisInds, s2t),...
%             initvx(thisInds, s2t), initvy(thisInds, s2t), height, width);
%         regWeights = params.lambda * (nbReg - repmat(directReg, ...
%             1, size(nb,2)));
%         cycleWeights = cycleWeights - regWeights;
%         cycleWeights(cycleWeights < 1) = 0;
%         weightsAll = nnW(thisInds, :) .* cycleWeights;
%         weightsAll = bsxfun(@times, weightsAll, 1./sum(weightsAll, 2));
%         pairvx(thisInds, s2t) = sum(nbvx .* weightsAll, 2);
%         pairvy(thisInds, s2t) = sum(nbvy .* weightsAll, 2);
%     end
%     if src > 1
%         % FIXME: get number of backspace needed automatically
%         fprintf(repmat('\b',1,7));
%     end
% end
% pairvx = round(pairvx);
% pairvy = round(pairvy);