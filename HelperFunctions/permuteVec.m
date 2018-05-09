function [ permutedVector ] = permuteVec( vecLength, chunkQuant )
% Permutes vector indexes based on number of windows (chunks) input
%   permutedVector: rearranged index vector (size [1xN])
%   vecLength:  length of desired permutation vector (size [1xN])
%   chunkQuant: number of windows over which to rearrange (scalar, max = N)

chunkSize = vecLength/chunkQuant;
if mod(chunkSize,1) ~= 0
    %extraTime = (chunkSize - floor(chunkSize))*chunkQuant;
    chunkSize = floor(chunkSize);
%else
    %extraTime = 0;
end

chunks = zeros(chunkQuant, chunkSize);
timeIndex = 1:vecLength;

for i=1:chunkQuant
    chunks(i,:) = timeIndex((i-1)*chunkSize+1 : i*chunkSize);
end

permutedVector = timeIndex; % leftover values are at end of time series data
randomChunk = randperm(chunkQuant);

for i=1:chunkQuant
    permutedVector((i-1)*chunkSize+1 : i*chunkSize) = chunks(randomChunk(i),:);
end

% use extraTime if needed?

end


