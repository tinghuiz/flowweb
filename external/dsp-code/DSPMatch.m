function [vx,vy,match_cost] = DSPMatch(sift1, sift2)
%% DSPMatch: Deformable spatial pyramid matching
% Match features (e.g., SIFTs) between two images
%
% Input:
% sift1: features from image1: h1 x w1 x n matrix
% (h1: image1's height, w1: image1' width, n: feature dimension)
% sift1: features from image2: h2 x w2 x n matrix
% (h2: image2's height, w2: image2' width, n: feature dimension)
%
% Output:
% vx (h1 x w1 matrix): x-translation of pixels in image1 to their corresponding pixels in image2 (x2 = x1 + vx)
% vy (h1 x w1 matrix): y-translation of pixels in image1 to their corresponding pixels in image2 (y2 = y1 + vy)
% match_cost (h1 x w1 matrix) = |sift1(y1,x1) - sift2(y2,x2)|
%
% Citation:
% Deformable Spatial Pyramid Matching for Fast Dense
% Correspondences, Jaechul Kim, Ce Liu, Fei Sha, and Kristen Grauman.
% CVPR 2013
%
% Author:
% Jaechul Kim, jaechul@cs.utexas.edu
 
 
[h1, w1, ~] = size(sift1);
[h2, w2, ~] = size(sift2);

%% Grid-layer
% parameters
n_init_target_feat = 2500;
init_grid_size = floor(sqrt(h1*w1/n_init_target_feat) + 0.5);
params.grid_size = init_grid_size;
params.n_tree_level = 3;
sx = max(w1, w2);
sy = max(h1, h2);
params.search_radius = [ceil(sx/2), ceil(sy/2)]; 
params.truncate_const = 500;

% data-term for each grid cell
[node_cost, node_trans, aux_info] = SPNodeMatchMex(sift1, sift2, params);
for i = 1 : params.n_tree_level
    offset = (4^(i-1)-1)/3;
    s_idx = offset + 1;
    e_idx = offset + 4^(i-1);
    x = node_cost(:, s_idx:e_idx);
    node_cost(:, s_idx:e_idx) = node_cost(:, s_idx:e_idx)./mean(x(:));
end

% BP
bp_params.deformation_coeff = 0.5*0.01;
bp_params.truncate_const = 0.25;
bp_params.max_iter = 50;
bp_params.n_tree_level = params.n_tree_level;
bp_params.n_xstate = aux_info(1);
bp_params.n_ystate = aux_info(2);
opt_state = SPBPMex(node_cost, bp_params);

%% Pixel-layer
% initialize from grid-layer solution
n_leaf = 4^(params.n_tree_level-1);
node_disp = node_trans(:, opt_state+1);
leaf_node_disp = node_disp(:, end-n_leaf+1:end);
n_bin = 2^(params.n_tree_level-1);
x = linspace(1, w1, n_bin+1);
y = linspace(1, h1, n_bin+1);
x = floor(x);
y = floor(y);
pixel_disparity_x = zeros(h1,w1, 'int32');
pixel_disparity_y = zeros(h1,w1, 'int32');
for i = 1 : n_bin
    for j = 1 : n_bin
        node_idx = j + n_bin*(i-1);
        lx = x(i);
        rx = x(i+1);
        ty = y(j);
        dy = y(j+1);
        pixel_disparity_x(ty:dy,lx:rx) = leaf_node_disp(1, node_idx);
        pixel_disparity_y(ty:dy,lx:rx) = leaf_node_disp(2, node_idx);
    end
end
in_disp = zeros(2, w1*h1, 'int32');
in_disp(1,:) = pixel_disparity_x(:)';
in_disp(2,:) = pixel_disparity_y(:)';

% pixel match: two-level coarse-to-fine search for speed-up
pxm_params.truncate_const = params.truncate_const;
pxm_params.search_grid_size = 4;
sx = 20;
sy = 20;
pxm_params.search_radius = floor([sx, sy]);
pxm_params.deform_coeff = 2.5;
out_disp = PixelMatchMex(sift1, sift2, in_disp, pxm_params);

pxm_params.search_grid_size = 1;
sx = 4;
sy = 4;
pxm_params.search_radius = floor([sx, sy]);
pxm_params.deform_coeff = 2.5;
[out_disp, match_cost] = PixelMatchMex(sift1, sift2, out_disp, pxm_params);

% output
vx = reshape(out_disp(1,:), [h1, w1]);
vy = reshape(out_disp(2,:), [h1, w1]);
vx = double(vx);
vy = double(vy);
match_cost = reshape(match_cost, [h1, w1]);

