function [ outputData ] = homeDetrend( inputData )
%HomeDetrend simply subtracts the ensemble mean and divides by the ensemble
%standard deviation for all realizations of inputData (effectively making
%it weakly stationary)
%   Input dims:  Channel/layer x time x trial

ensembleMean = repmat(mean(inputData,3),1,1,size(inputData,3));
ensembleStdev = repmat(std(inputData,0,3),1,1,size(inputData,3));
outputData = (inputData-ensembleMean)./ensembleStdev;

end

