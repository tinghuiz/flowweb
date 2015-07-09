function [accuracy, seg] = EvaluateCT101Match(anno1, anno2, vx, vy)

[h1 w1] = size(anno1);
[h2 w2] = size(anno2);

[x1 y1] = meshgrid(1:w1, 1:h1);
x2 = x1+ vx;
y2 = y1+ vy;

% localization
[obj_x1 obj_y1] = ObjPtr(anno1);
[obj_x2 obj_y2] = ObjPtr(anno2);

tf = x2 >= 1 & x2 <= w2 & y2 >= 1 & y2 <= h2;
ptr_ind1 = sub2ind(size(anno1), y1(tf), x1(tf));
ptr_ind2 = sub2ind(size(anno2), y2(tf), x2(tf));
loc_err_map = inf(size(anno1));
loc_err_map(ptr_ind1) =...
    abs(obj_x1(ptr_ind1) - obj_x2(ptr_ind2)) + abs(obj_y1(ptr_ind1) - obj_y2(ptr_ind2));

% localization evaluation
[seg, in_bound] = TransferAnnotation(anno2, vx,vy);   
% true_match = seg == anno1 & anno1 == 1 & in_bound;
% loc_err.correct_fg = mean2(loc_err_map(true_match));
% fg_match = seg == 1 & in_bound;
% loc_err.fg = mean2(loc_err_map(fg_match));
% loc_err.all = mean2(loc_err_map(in_bound));
loc_err = mean2(loc_err_map(in_bound));

% label transfer evluation
mean_acc = mean2(seg == anno1);
accuracy.mean = mean_acc;
% fg = sum(seg(:) == anno1(:) & anno1(:) == 1)./sum(anno1(:));
% bg = sum(seg(:) == anno1(:) & anno1(:) == 0)./sum(anno1(:)==0);
%accuracy.fg = fg;
%accuracy.bg = bg;
i = seg == 1 & anno1 == 1;
u = seg == 1 | anno1 == 1;
accuracy.iou = sum(i(:))/sum(u(:));
accuracy.loc_err = loc_err;

end

function [x1 y1] = ObjPtr(anno)
[y1 x1] = find(anno);
lx1 = min(x1);
rx1 = max(x1);
ty1 = min(y1);
dy1 = max(y1);
w1 = rx1-lx1 + 1;
h1 = dy1-ty1 + 1;
[x y] = meshgrid(1:size(anno,2), 1:size(anno,1));
x1 = (x - lx1)./w1;
y1 = (y - ty1)./h1;
end
