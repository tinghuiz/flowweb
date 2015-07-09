function D = fdist(vx1, vy1, vx2, vy2, height, width)
% Per-flow distance between two flow fields
vx1 = double(vx1)/width;
vy1 = double(vy1)/height;
vx2 = double(vx2)/width;
vy2 = double(vy2)/height;
f1 = [vx1(:), vy1(:)];
f2 = [vx2(:), vy2(:)];
D = sqrt(sum((f1 - f2).^2, 2));
D = reshape(D, size(vx1));
