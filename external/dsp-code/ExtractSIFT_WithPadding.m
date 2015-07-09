function [sift1, mag1] = ExtractSIFT_WithPadding(im1, pca_basis, sift_size)

if ( nargin <= 2)
    sift_size = 4;
end

if  (ndims(im1) == 3 )
    im1 = rgb2gray(im1);
end
im1 = single(im1);


% pad
pad_size = ceil(1 + 3/2 * sift_size);
[h0,w0,~] = size(im1);
im1 = padarray(im1, [pad_size, pad_size], 'replicate', 'both');
pad_lx = pad_size + 1;
pad_ty = pad_size + 1;

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
ty = min(f1(2,:));
blx = pad_lx - lx + 1;
bty = pad_ty - ty + 1;
    
% cut -- back to original image
sift1 = sift1(bty : bty + h0-1, blx : blx + w0-1, :);
mag1 = mag1(bty : bty + h0-1, blx : blx + w0-1);

