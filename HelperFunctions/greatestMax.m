function [ maxOfMaxes ] = greatestMax( matrix )
% greatestMax finds the maximum over all dimensions of a matrix 
%   matrix = your input matrix; can be as many dimensions as desired
%   maxOfMaxes = the highest single value in matrix

% for some odd reason, if you make this 1 instead of 2 it gets stuck...
% this works with the additional line at the bottom, though I'm not sure 
% why I need to do it this way? Probably just a minor bug on my end.

while length(size(matrix)) > 2
    matrix = squeeze(max(matrix));
end

maxOfMaxes = max(max(matrix));

