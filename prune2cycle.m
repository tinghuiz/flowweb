function pmask = prune2cycle(pairvx, pairvy, gx, gy, epsilon)
% 2-cycle check for the input FlowWeb

N = sqrt(size(pairvx, 2));
pmask = false(size(pairvx));
for src = 1 : N
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        s2t = src + (tgt - 1) * N;
        t2s = tgt + (src - 1) * N;
        pmask(:, s2t) = pruneOnePair(pairvx(:,s2t), pairvy(:,s2t), ...
            pairvx(:,t2s), pairvy(:,t2s), gx, gy, epsilon);
    end
end

function pmask = pruneOnePair(vx1, vy1, vx2, vy2, gx, gy, epsilon)

[H, W] = size(gx);
pmask = false(H*W, 1);
px12 = round(gx(:) + vx1(:));
py12 = round(gy(:) + vy1(:));
invalid = px12 > W | px12 < 1 | py12 > H | py12 < 1;
px12(invalid) = 1;
py12(invalid) = 1;
idx = py12(:) + (px12(:) - 1) * H;
px121 = px12(:) + vx2(idx);
py121 = py12(:) + vy2(idx);
% FIXME: actual euc distance?
xdiff = abs(px121 - gx(:));
ydiff = abs(py121 - gy(:));
pmask(invalid | xdiff > epsilon | ydiff > epsilon) = true;
