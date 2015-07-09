function ims = loadImageSet(imdir)

imfiles = dir([imdir, '*.jpg']);
N = length(imfiles);

ims = cell(1, N);
for i = 1 : N
    ims{i} = im2double(imread([imdir, imfiles(i).name]));
    if size(ims{i},3) == 1
        ims{i} = repmat(ims{i}, [1 1 3]);
    end
end

