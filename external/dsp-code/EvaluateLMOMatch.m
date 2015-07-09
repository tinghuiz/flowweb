function [accuracy, seg] = EvaluateLMOMatch(anno1, anno2, vx, vy)

[seg, in_bound] = TransferAnnotation(anno2, vx,vy);   

% fill out-bound label
[d l] = bwdist(in_bound);
seg(~in_bound) = seg(l(~in_bound));

% evaluation
label1 = setdiff(unique(anno1(:)), 0);
label2 = setdiff(unique(anno2(:)), 0);
common_label = intersect(label1,label2);
tf = false(size(anno1));
for i = 1 : numel(common_label)
    tf(anno1 == common_label(i)) = true;
end
mean_acc = mean2(seg(tf) == anno1(tf));
accuracy.mean = mean_acc;
accuracy.n_pixel = sum(tf(:));
accuracy.n_hit = sum(seg(tf) == anno1(tf));
