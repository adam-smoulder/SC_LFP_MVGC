%% Channel isolation script
% takes rawData -> rawDataGoodChan for the data variable
% only pulls out the selected channels for analysis
% -load data variable before running

% get channel index order rearranged:
channelOrder = data(1).channelOrder;
channelIndex = zeros(size(channelsToUse));
for i=1:numel(channelsToUse)
    channelIndex(ceil(i/2),2-mod(i,2)) = ...
        find(channelsToUse(ceil(i/2),2-mod(i,2)) == channelOrder);
end

flatChanIndex = reshape(channelIndex',[1 6]); % flatten for quick check

% only select desired channels
for i=1:length(data)
       data(i).rawDataGoodChan = ...
           double(data(i).rawData(flatChanIndex,:));
end