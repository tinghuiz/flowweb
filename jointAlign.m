function [pairvx, pairvy] = jointAlign(pairvx, pairvy, params)


N = size(pairvx,1);
H = size(pairvx{1,2}, 1);
W = size(pairvx{1,2}, 2);
params.height = H;
params.width = W;
params.numImage = N;
epsilon = params.cycleThresh * max([H, W]);
[pairvx, pairvy] = flowCellToMat(pairvx, pairvy);
[gx, gy] = meshgrid(1:W, 1:H);
initvx = pairvx;
initvy = pairvy;

% Pre-compute spatial weights for the intra-image phase
spatialWeights = slmetric_pw([gx(:), gy(:)]', [gx(:), gy(:)]', 'sqdist');
spatialSigma = params.filterSpatialSigma * max([H, W]);
spatialWeights = exp(-spatialWeights/spatialSigma);
spatialWeights = single(spatialWeights);
% FIXME: this thresholding necessary?
spatialWeights(spatialWeights < 1e-3) = 0;

iter = 1;
while 1
    %% Inter-image phase:
    fprintf('Iter %d, Inter-image Phase:\n', iter);
    % 2-cycle check
    fprintf('Checking 2-cycles...\n');
    pruneMask = prune2cycle(pairvx, pairvy, gx, gy, epsilon);
    % Get 3-cycle set
    fprintf('Computing 3-cycle set...\n');
    cycleSet = eval3cycle(pairvx, pairvy, gx, gy, pruneMask, epsilon);
    % Prune ones that did not pass the 2-cycle check
    cycleSet = cycleSet & repmat(~pruneMask, [1, 1, N]);
    % Compute objective
    obj = objective(cycleSet, pairvx, pairvy, initvx, initvy, H, W, params.lambda);
    % Check convergence
    if iter > 1
        thisObj = obj;
        improvePercent = (thisObj - prevObj)/prevObj;
        fprintf('Prev obj. = %.3f, this obj. = %.3f, relative improvement: %.3f\n',...
            prevObj, thisObj, improvePercent);
        if improvePercent < params.convergeThresh
            break;
        end
    end
    prevObj = obj;
    % Compute priority
    fprintf('Computing flow update priority...\n');
    [priorScores, transits] = priority(pairvx, pairvy, gx, gy, cycleSet,...
        initvx, initvy, params.lambda);
    % Update flows
    fprintf('Updating ....\n');
    [pairvx, pairvy] = updateflows(pairvx, pairvy, priorScores,...
        transits, params);

    %% Intra-image phase:
    fprintf('Iter %d, Intra-image Phase:\n', iter);
    % 2-cycle check
    fprintf('Checking 2-cycles...\n');
    pruneMask = prune2cycle(pairvx, pairvy, gx, gy, epsilon);
    % Get 3-cycle set
    fprintf('Computing 3-cycle set...\n');
    cycleSet = eval3cycle(pairvx, pairvy, gx, gy, pruneMask, epsilon);
    % Prune cycle_set with prune_mask
    cycleSet = cycleSet .* repmat(~pruneMask, [1, 1, N]);
    fprintf('Consistency-weighted filtering...\n');
    [pairvx, pairvy] = cyclefilter(pairvx, pairvy, cycleSet, ...
        spatialWeights, params, true);
    if params.saveEveryIter
        pvx = pairvx;
        pvy = pairvy;
        pairvx = int16(round(pairvx));
        pairvy = int16(round(pairvy));
        [pairvx, pairvy] = flowMatToCell(pairvx, pairvy, H);
        fname = sprintf('%s/pairflows_iter_%.2d.mat', ...
            params.resDir, iter);
        save(fname, 'pairvx', 'pairvy');
        pairvx = pvx;
        pairvy = pvy;
    end
    iter = iter + 1;
    if iter > params.maxIter
        break;
    end
end