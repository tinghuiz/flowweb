function [sift1, bbox1, mag1] = ExtractSIFT(im1, pca_basis, sift_size)

if ( nargin <= 2)
    sift_size = 4;
end

if  (ndims(im1) == 3 )
    im1 = rgb2gray(im1);
end
im1 = single(im1);

% Extract SIFT
[f1, sift1] = vl_dsift(im1, 'step', 1, 'size', sift_size, 'norm', 'FloatDescriptors', 'fast');

if ( ~isempty(pca_basis) ) % dimensionality reduction
    sift1 = pca_basis'*sift1;
end

% Formatting
x1 = unique(f1(1,:));
y1 = unique(f1(2,:));
width1 = numel(x1);
height1 = numel(y1);
sift1 = reshape(sift1', [height1, width1, size(sift1,1)]);
mag1 = reshape(f1(3,:), [height1, width1]);

lx = min(f1(1,:));
rx = max(f1(1,:));
ty = min(f1(2,:));
dy = max(f1(2,:));
bbox1 = [lx,rx,ty,dy]';
