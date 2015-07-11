function res = evaluate(pairvx, pairvy, annos, params)

switch params.metric
    case 'keypoint'
        res = pck(pairvx, pairvy, annos, params);
    case 'part'
        res = partIOU(pairvx, pairvy, annos, params);
end

function [pckMean, pckAll] = pck(pairvx, pairvy, annos, params)
N = params.numImage;
W = params.width;
H = params.height;
thresh = params.alpha * max([W, H]);
pcnt = 0;
pckAll = zeros(N, N);
for src = 1 : N
    for tgt =1 : N
        if src == tgt
            continue;
        end
        srcKeypts = annos{src}.keypts;
        tgtKeypts = annos{tgt}.keypts;
        if isempty(srcKeypts) || isempty(tgtKeypts)
            continue;
        end
        valid_mask = annos{src}.keyptStatus & annos{tgt}.keyptStatus & ...
            srcKeypts(1,:) < W & srcKeypts(1,:) > 1 & ...
            srcKeypts(2,:) < H & srcKeypts(2,:) > 1 & ...
            tgtKeypts(1,:) < W & tgtKeypts(1,:) > 1 & ...
            tgtKeypts(2,:) < H & tgtKeypts(2,:) > 1;
        if sum(valid_mask) == 0
            continue;
        end
        pcnt = pcnt + 1;
        srcKeypts = srcKeypts(:, valid_mask);
        tgtKeypts = tgtKeypts(:, valid_mask);
        srcKeyInds = round(srcKeypts(2,:)) + (round(srcKeypts(1,:)) - 1) * H;
        tpx = srcKeypts(1,:) + double(pairvx{src, tgt}(srcKeyInds));
        tpy = srcKeypts(2,:) + double(pairvy{src, tgt}(srcKeyInds));
        error = sqrt((tpx - tgtKeypts(1,:)).^2 + (tpy - tgtKeypts(2,:)).^2);
        correct = error <= thresh;
        pckAll(src,tgt) = sum(correct)/numel(error);
    end
end
pckMean = sum(pckAll(:))/pcnt;

function [iouMean, iouAll] = partIOU(pairvx, pairvy, annos, params)
N = params.numImage;
W = params.width;
H = params.height;
iouAll = zeros(N, N);
for src = 1 : N
    fprintf('Eval parts: %d/%d\n', src, N);
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        unionParts = [];
        iouParts = [];
        for sp = 1 : length(annos{src}.parts)
            for tp = 1 : length(annos{tgt}.parts)
                if strcmp(annos{tgt}.parts(tp).part_name, ...
                        annos{src}.parts(sp).part_name)
                    tgtMask = annos{tgt}.parts(tp).mask;
                    srcMask = annos{src}.parts(sp).mask;
                    [warpTgtMask, ~] = warpImage(single(tgtMask), ...
                        single(pairvx{src, tgt}), ...
                        single(pairvy{src, tgt}));
                    warpTgtMask(isnan(warpTgtMask)) = 0;
                    union = sum(sum(warpTgtMask | srcMask));
                    inter = sum(sum(warpTgtMask & srcMask));
                    iouParts = [iouParts; inter/(union + eps)];
                    unionParts = [unionParts; union];
                end
            end
        end
        weightedIOU = iouParts .* (unionParts./(sum(unionParts) + eps));
        iouAll(src, tgt) = sum(weightedIOU);
    end
end
iouMean = sum(iouAll(:)) / (N^2 - N);