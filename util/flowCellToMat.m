function [mvx, mvy] = flowCellToMat(cvx, cvy)
% Convert the FlowWeb from cells to a matrix such that cvx{i,j}(:) is 
% stored in mvx(:, i+(j-1)*N)

N = size(cvx,1);
M = numel(cvx{1,2});
mvx = single(zeros(M, N^2));
mvy = single(zeros(M, N^2));
for src = 1 : N
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        s2t = src + (tgt - 1) * N;
        mvx(:, s2t) = single(cvx{src, tgt}(:));
        mvy(:, s2t) = single(cvy{src, tgt}(:));
    end
end