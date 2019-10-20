function Y = BuildY(img, blockSize)

% function Y = BuildY(img, blockSize)
% Input:
%	img - the image to decompose into patches
%	blockSize - the patch size (assumes square patches)
%
% Output:
%	Y - the patch matrix, patches are in rows
%
% Note that the image will be cropped to the nearest multiple of blockSize

if ~exist('blockSize', 'var') || isempty(blockSize)
    blockSize = 8;
end

% BUILDY Take an input image img and break it down into 8x8 patches and make
% a column vector out of each of these patches

img = img(1:floor(end/blockSize)*blockSize, 1:floor(end/blockSize)*blockSize, :);

Y = zeros(size(img, 1)/blockSize * size(img, 2)/blockSize, blockSize^2 * size(img, 3), class(img));

for i = 1:blockSize:size(img, 2)
    for c = 1:size(img, 3)
        t = img(:, i:i+blockSize-1, c)';
        Y((i-1)/blockSize*size(img, 1)/blockSize+1:(i-1)/blockSize*size(img, 1)/blockSize+size(t, 2)/blockSize, (c - 1)*blockSize^2+1:c*blockSize^2) = ...
            reshape(t, [blockSize^2 size(t, 2)/blockSize])';
    end
end

end