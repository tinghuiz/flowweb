function [cvx, cvy] = flowMatToCell(mvx, mvy, imgHeight)
% Convert the FlowWeb from a matrix to cells such that mvx(:, i+(j-1)*N)
% is stored in cvx{i,j}(:)

N = sqrt(size(mvx, 2));
cvx = cell(N, N);
cvy = cell(N, N);
for src = 1 : N
    for tgt = 1 : N
        if src == tgt
            continue;
        end
        s2t = src + (tgt - 1) * N;
        cvx{src,tgt} = reshape(mvx(:, s2t), imgHeight, []);
        cvy{src,tgt} = reshape(mvy(:, s2t), imgHeight, []);
    end
end