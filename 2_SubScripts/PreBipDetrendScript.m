%% Pre-Bip detrend
% runs homeDetrend on data before the bipolar subtraction occurs
% get data in one big matrix
tempTargData = zeros([size(data(1).targlfpmat) length(data)]);
tempSaccData = zeros([size(data(1).sacclfpmat) length(data)]);
for i=1:length(data)
    tempTargData(:,:,i) = data(i).targlfpmat;
    tempSaccData(:,:,i) = data(i).sacclfpmat;
end

% detrend
tempTargData = homeDetrend(tempTargData);
tempSaccData = homeDetrend(tempSaccData);

% reassign data
for i=1:length(data)
    data(i).targlfpmat = squeeze(tempTargData(:,:,i));
    data(i).sacclfpmat = squeeze(tempSaccData(:,:,i));
end

disp('Finished detrending')

if showDetrendSgrams
    figureCount = figureCount+1;
    figure(figureCount)
    SgramFromAlignedData
    hold on
    subplot(3,2,1)
    title(['Pre - detrend, trial ' num2str(trialNum)])
end