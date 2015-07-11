function [priorScores, transits] = priority(pairvx, pairvy, gx, gy, ... 
    cycleSet, initvx, initvy, lambda)
% Compute priority score for each flow, and record the transitive image
% that achieves the maximum

[H, W] = size(gx);
N = sqrt(size(pairvx, 2));
% Pre-allocate
hubPixelInds = int32(zeros(size(pairvx)));
transitCycleBound = zeros(H*W, N);
tvx = zeros(H*W, N);
tvy = zeros(H*W, N);
outMask = false(size(pairvx));
priorScores = single(zeros(size(pairvx)));
transits = single(zeros(size(pairvx)));
% Pre-compute pixel indices in the hub/transition image
for src = 1 : N
    for hub = 1 : N
        if src == hub
            continue;
        end
        s2h = src + (hub - 1) * N;
        hpx = gx(:) + pairvx(:, s2h);
        hpy = gy(:) + pairvy(:, s2h);
        outMask(:, s2h) = hpx > W | hpx < 1 | hpy > H | hpy < 1;
        hpx(outMask(:, s2h)) = 1;
        hpy(outMask(:, s2h)) = 1;
        hubPixelInds(:, s2h) = int32(hpy(:) + (hpx(:) - 1) * H);
    end
end

pcnt = 0;
for src = 1 : N
    % Assuming 2-cycle/symmetry holds, only need to compute half the pairs
    for tgt = src+1 : N
        pcnt = pcnt + 1;
        if pcnt > 1
            fprintf('%.6d/%.6d', pcnt, (N^2 - N)/2);
        end
        transitCycleBound(:) = 0; 
        for hub = 1 : N
            if src == hub || tgt == hub
                continue;
            end
            s2h = src + (hub - 1) * N;
            h2t = hub + (tgt - 1) * N;
            cycle_sh  = squeeze(cycleSet(:, s2h, :));
            cycle_ht  = squeeze(cycleSet(:, h2t, :));
            cycle_sht = cycle_sh & cycle_ht(hubPixelInds(:, s2h), :);
            transitCycleBound(:, hub) = sum(cycle_sht, 2)/(N-2);
            % Compute transitive flows
            tvx(:, hub) = pairvx(:, s2h) + pairvx(hubPixelInds(:, s2h), h2t);
            tvy(:, hub) = pairvy(:, s2h) + pairvy(hubPixelInds(:, s2h), h2t);
            % FIXME: is the following necessary? for out of range flows,
            % the cycle_sh is already zero?
%             hub_cycle_bound(out_range(:, pid_sh), hub) = 0;
        end
        s2t = src + (tgt - 1) * N;
        directCycle = sum(squeeze(cycleSet(:, s2t, :)), 2)/(N-2);
        transitReg = fdist(tvx, tvy, repmat(initvx(:, s2t), 1, N), ...
            repmat(initvy(:, s2t), 1, N), H, W);
        directReg = fdist(pairvx(:, s2t), pairvy(:,s2t), initvx(:, s2t),...
            initvy(:, s2t), H, W);
        hubScores = transitCycleBound - repmat(directCycle(:), [1, N]) - ...
            lambda * (transitReg - repmat(directReg, 1, N));
        [priorScores(:, s2t), maxTransitInds] = max(hubScores, [], 2);
        transits(:, s2t) = maxTransitInds;
        if pcnt > 1
            % FIXME: get number of backspace needed automatically
            fprintf(repmat('\b',1,13));
        end
    end
end