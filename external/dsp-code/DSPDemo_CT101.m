% input
% im1, im2 : input images to match
% anno1, anno2: ground truth object annotation for evaluation
% pca_basis: pca basis for dimensionality reduction of sift

load ct101_example_data.mat im1 im2 anno1 anno2
load ct101_pca_basis.mat pca_basis 
% pca_basis = []; % if you want to use the sift of original dimension
sift_size = 4; % sift patch size 16 x 16 pixels (i.e., patch_size = 4*sift_size)

% extract SIFT
[sift1, bbox1] = ExtractSIFT(im1, pca_basis, sift_size);
[sift2, bbox2] = ExtractSIFT(im2, pca_basis, sift_size);
im1 = im1(bbox1(3):bbox1(4), bbox1(1):bbox1(2), :);
im2 = im2(bbox2(3):bbox2(4), bbox2(1):bbox2(2), :);
anno1 = anno1(bbox1(3):bbox1(4), bbox1(1):bbox1(2), :);
anno2 = anno2(bbox2(3):bbox2(4), bbox2(1):bbox2(2), :);

% Match
tic; 
[vx,vy] = DSPMatch(sift1, sift2); 
t_match = toc;

% Evaluation
[acc, seg] = EvaluateCT101Match(anno1, anno2, vx, vy);
acc.time = t_match;

% Warping
warp21=warpImage(im2double(im2),vx,vy); % im2 --> im1

disp('----------------------------')
disp('accuracy')
disp(acc)

figure,
subplot(2,2,1);
imshow(im1);
title('image1');
subplot(2,2,3);
imshow(im2);
title('image2');
subplot(2,2,2);
imshow(warp21);
title('warp 2-->1');
subplot(2,2,4);
imshow(seg);
title('label transfer 2-->1');