function [nbInds, nbWeights] = setupSpatialFilter(params)
% Return nearby pixel indices and weights for use in the intra-image phase

W = params.width;
H = params.height;
[gx, gy] = meshgrid(1:W, 1:H);
dist = pdist2([gx(:), gy(:)], [gx(:), gy(:)]);
sigma = params.filterSpatialSigma * max([H, W]);
weights = exp(-dist/sigma);
% FIXME: replace with bsxfun for speed up
weights = weights ./ repmat(sum(weights, 2), [1, H*W]);
weights = single(weights);
% Only consider pixels within certain radius
weights(weights < params.spatialWeightThresh) = 0;
% Renormalize again
weights = weights ./ repmat(sum(weights, 2), [1, H*W]);
nbInds = cell(H*W, 1);
nbWeights = cell(H*W,1);
maxlen = 0;
for i = 1 : length(nbInds)
    nbInds{i}  = find(weights(i,:));
    nbWeights{i} = weights(i, nbInds{i});
    if numel(nbInds{i}) > maxlen
        maxlen = numel(nbInds{i});
    end
end
nbInds = cell2mat(cellfun(@(x) cat(2, x, ones(1,maxlen-length(x))),...
    nbInds,'UniformOutput',false));
nbWeights = cell2mat(cellfun(@(x) cat(2, x, zeros(1,maxlen-length(x))),...
    nbWeights,'UniformOutput',false));