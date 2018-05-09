% used in preprocessing to incorporate data of all channels

for trial = 1:nTrials % for each trial in the data
    disp(['Bipolar Subtraction for Trial ' num2str(trial)])
    
    % load data
    targLfp = data(trial).targlfpmat(channelsToUse, :);
    saccLfp = data(trial).sacclfpmat(channelsToUse, :);
    
    bipTargLfp = zeros(size(targLfp));
    bipSaccLfp = zeros(size(saccLfp));
    
    % superficial
    bipTargLfp(1,:) = ((targLfp(supChans(1),:)-targLfp(supChans(3),:)) + ...
        (targLfp(supChans(4),:)-targLfp(supChans(4),:)))/2;
    bipSaccLfp(1,:) = ((saccLfp(supChans(1),:)-saccLfp(supChans(3),:)) + ...
        (saccLfp(supChans(4),:)-saccLfp(supChans(4),:)))/2;
    
    % intermediate
    bipTargLfp(2,:) = ((targLfp(intChans(1),:)-targLfp(intChans(3),:)) + ...
        (targLfp(intChans(4),:)-targLfp(intChans(4),:)))/2;
    bipSaccLfp(2,:) = ((saccLfp(intChans(1),:)-saccLfp(intChans(3),:)) + ...
        (saccLfp(intChans(4),:)-saccLfp(intChans(4),:)))/2;
    
    % deep
    bipTargLfp(3,:) = ((targLfp(deepChans(1),:)-targLfp(deepChans(3),:)) + ...
        (targLfp(deepChans(4),:)-targLfp(deepChans(4),:)))/2;
    bipSaccLfp(3,:) = ((saccLfp(deepChans(1),:)-saccLfp(deepChans(3),:)) + ...
        (saccLfp(deepChans(4),:)-saccLfp(deepChans(4),:)))/2;
    
    % save results to data
    data(trial).targlfpbip = bipTargLfp;
    data(trial).sacclfpbip = bipSaccLfp;
    
end