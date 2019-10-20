function img = YToImage(Y, dim, blockSz)

% function img = YToImage(Y, dim, blockSz)
% Input:
%   Y - the patch matrix (patches are in rows)
%   dim - the size of the original image (will be rounded automatically
%	      to the nearest whole multiple of patch size)
%   blockSz - the block size
%         If patch size is a multiple of 3, will be assumed to be a color image
%         Default is sqrt(size(Y, 2)) for grayscale (patch size does not divide evenly by 3)
%         and sqrt(size(Y, 2) / 3) for color.
%         set dim(3) = 3
%		  
% Output:
%	img - the reconstructed image
%

blockCnt = size(Y, 1);

if ~exist('blockSz', 'var') || isempty(blockSz)
    % Asssume a color image if patch size divides by 3
    if mod(size(Y, 2), 3) == 0
        blockSz = sqrt(size(Y, 2) / 3);
    else
        blockSz = sqrt(size(Y, 2));
    end
end

if blockSz ~= round(blockSz)
	error('blockSz has to be an integer')
end

columnBlocks = floor(dim(1) / blockSz);
dimY = columnBlocks * blockSz;

rowBlocks = blockCnt / columnBlocks;

if rowBlocks ~= round(rowBlocks)
    error('Image dimension doesn`t match');
end

% For compatibility with previous definition
if length(dim) < 3
    dim(3) = 1;
end

dimX = rowBlocks * blockSz;
img = zeros(dimY, dimX, dim(3));

for i = 1:blockSz:dimX
    for c = 1:dim(3)
		img(:, i:i+blockSz-1, c) = reshape( ...
			Y(columnBlocks * (i - 1) / blockSz + 1:columnBlocks * ((i - 1) / blockSz + 1), ...
				(c - 1) * blockSz^2 + 1:c * blockSz^2)', ...
			[blockSz dimY])';
    end
end

end