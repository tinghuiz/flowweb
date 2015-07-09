function [pairvx, pairvy] = updateflows(pairvx, pairvy, priorScores, ...
    transits, params)

%% Sort the priority queue
N = sqrt(size(pairvx, 2));
M = size(pairvx, 1);
H = params.imgHeight;
numKeep = round(M * params.proposalRate);
updateQueue = zeros((N^2-N)/2 * numKeep, 5);

first = 1;
for src = 1 : N
    for tgt = src+1 : N
        last = first + numKeep -1;
        s2t = src + (tgt - 1) * N;
        % Only keep top 'numKeep' flows for consideration
        [topScores, topInds] = maxk(priorScores(:, s2t), numKeep);
        updateQueue(first:last, 1) = src;
        updateQueue(first:last, 2) = tgt;
        updateQueue(first:last, 3) = transits(topInds, s2t);
        updateQueue(first:last, 4) = topInds;
        updateQueue(first:last, end) = topScores;
        first = first + numKeep;
    end
end
updateQueue = updateQueue(1:last, :);
posInds = updateQueue(:,end) > 0;
updateQueue = updateQueue(posInds, :);
% FIXME: below necessary?
updateQueue(:,end) = updateQueue(:,end) + 1e-6 * randn(size(updateQueue,1), 1);
numUpdate = round(params.updateRate * size(updateQueue,1));
[~, topInds] = maxk(updateQueue(:,end), numUpdate);
updateQueue = updateQueue(topInds, :);

%% Actual update
src = updateQueue(:, 1);
tgt = updateQueue(:, 2);
hub = updateQueue(:, 3);
srcPixelInds = updateQueue(:, 4);
spy = rem(srcPixelInds - 1, H) + 1;
spx = (srcPixelInds - spy)/H + 1;
s2h = src + (hub - 1) * N;
h2t = hub + (tgt - 1) * N;
s2t = src + (tgt - 1) * N;
hpx = spx + pairvx((s2h - 1) * M + srcPixelInds);
hpy = spy + pairvy((s2h - 1) * M + srcPixelInds);
hubPixelInds = hpy + (hpx - 1) * H;
tpx = hpx + pairvx((h2t - 1) * M + hubPixelInds);
tpy = hpy + pairvy((h2t - 1) * M + hubPixelInds);

% Update src to tgt with transitive flows
pairvx((s2t - 1) * M + srcPixelInds) = tpx - spx;
pairvy((s2t - 1) * M + srcPixelInds) = tpy - spy;
% Update tgt to src
t2s = tgt + (src - 1) * N;
tgtPixelInds = tpy + (tpx - 1) * H;
% FIXME: matlab not happy without the 'unique' pruning
t2sUpdateInds = (t2s - 1) * M + tgtPixelInds;
[t2sUpdateInds, uniqueInds] = unique(t2sUpdateInds);
% FIXME: see if removing vx_t2s and vy_t2s possible
vx_t2s = spx - tpx;
vy_t2s = spy - tpy;
pairvx(t2sUpdateInds) = vx_t2s(uniqueInds);
pairvy(t2sUpdateInds) = vy_t2s(uniqueInds);
