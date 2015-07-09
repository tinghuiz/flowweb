function obj = objective(cycleSet, thisvx, thisvy, initvx, initvy, height, width, lambda)

N = sqrt(size(thisvx, 2));
M = size(thisvx,1);
AFCC = sum(sum(single(sum(cycleSet, 3))/(N-2)))/(N*(N-1)*M);
R = sum(fdist(thisvx(:), thisvy(:), initvx(:), initvy(:), height, width))/(N*(N-1));
obj = AFCC - lambda * R;