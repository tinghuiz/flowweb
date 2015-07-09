function cycleSet = eval3cycle(pairvx, pairvy, gx, gy, pmask, epsilon)

[H, W] = size(gx);
N = sqrt(size(pairvx, 2));
cycleSet = false(H*W, N^2, N);

% Make flows that fail 2-cycle check go out of bounds
% FIXME: nicer solution?
pairvx(pmask) = W + 9999;
pairvy(pmask) = H + 9999;

for src = 1 : N
    if src > 1
        fprintf([repmat('\b',1,7) '%.3d/%.3d'], src, N);
    end
    for hub = 1 : N
        if src == hub
            continue;
        end
        s2h = src + (hub - 1) * N;
        hpx = gx(:) + pairvx(:, s2h);
        hpy = gy(:) + pairvy(:, s2h);
        out = hpx > W | hpx < 1 | hpy > H | hpy < 1;
        hpx(out) = 1;
        hpy(out) = 1;
        hubPixelInds = hpy(:) + (hpx(:) - 1) * H;
        for tgt = 1 : N
            if src == tgt || tgt == hub
                continue;
            end
            h2t = hub + (tgt - 1) * N;
            s2t = src + (tgt - 1) * N;
            % Transitive flows
            tvx = pairvx(:, s2h) + pairvx(hubPixelInds, h2t);
            tvy = pairvy(:, s2h) + pairvy(hubPixelInds, h2t);
            % FIXME: actual euc dist?
            xdiff = abs(tvx - pairvx(:, s2t));
            ydiff = abs(tvy - pairvy(:, s2t));
            
            cycleSet(:, s2t, hub) = xdiff <= epsilon & ...
                ydiff <= epsilon & ~out;
            % Note: Trying to speed up by setting cycle_set{sid,tid} =
            % cycle_set{tid,sid} is not correct, as the onmask in the source
            % image is different than that of the target image (needs to
            % calculate the correspondences, which makes the speed-up useless)
        end
    end
end
fprintf('\n');