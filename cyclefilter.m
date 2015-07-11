function [vx, vy] = cyclefilter(vx, vy, cycle_set, ...
    spatial_weights, params, is_gpu)

N = sqrt(size(vx,2));
cycle_map = int16(sum(cycle_set, 3));
numfilter = round(params.filterRate * params.width * params.height);

if is_gpu
    cycle_map = gpuArray(cycle_map);
    spatial_weights = gpuArray(spatial_weights);
    vx = gpuArray(vx);
    vy = gpuArray(vy);
else
    vx = gather(vx);
    vy = gather(vy);
    spatial_weights = gather(spatial_weights);
end

for src = 1 : N
    fprintf('filtering: %d/%d\n', src, N);
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        pid_st = src + (tgt - 1) * N;
        [~, fidx] = mink(cycle_map(:, pid_st), numfilter);
        
        cycle_diff = bsxfun(@minus, cycle_map(:, pid_st)', ...
           cycle_map(fidx, pid_st));
        cycle_weights = exp(single(cycle_diff)/params.filterCycleSigma/(N-2));
        cycle_weights(cycle_weights < 1) = 0;
        weights_all = spatial_weights(fidx, :) .* cycle_weights;
        weights_all = bsxfun(@times, weights_all, 1./sum(weights_all, 2));
        vx(fidx, pid_st) = sum(bsxfun(@times, weights_all, vx(:, pid_st)'),2);
        vy(fidx, pid_st) = sum(bsxfun(@times, weights_all, vy(:, pid_st)'),2);
    end
end
vx = round(vx);
vy = round(vy);
vx = gather(vx);
vy = gather(vy);

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