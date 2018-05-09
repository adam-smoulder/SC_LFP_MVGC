% used to make sure channels for neural data show activity
% load raw data beforehand

%% 1) sorting channels into order
sortedDataForCheck = zeros(3, 16, 1001, length(data)); % dims: cue, channel, time, trial

for trial = 1:length(data) % adjust for channel order
    trialData = data(trial);
    for i=1:length(trialData.channelOrder)
        sortedDataForCheck(1, trialData.channelOrder(i),:,trial) = trialData.targspikemat(i,:);
        sortedDataForCheck(2, trialData.channelOrder(i),:,trial) = trialData.gospikemat(i,:);
        sortedDataForCheck(3, trialData.channelOrder(i),:,trial) = trialData.saccspikemat(i,:);
    end
end

%% 2) plot data and show target value range

timeVec = [-400:1:600 ; -400:1:600 ; -600:1:400];
cue = {'targ', 'go', 'sacc'};
channelMax = NaN(3,16);

for i = 1:size(sortedDataForCheck,2) % for each channel
    figure()
    for j = 1:3 % for each cue (targ, go, sacc)
        % plotting
        subplot(3,1,j)
        hold on
        y = squeeze(mean(sortedDataForCheck(j,i,:,:),4));
        cuePoint = -timeVec(j,1);
        plot(timeVec(j,:), y)
        
        % evaluate average value near cue
        channelMax(j,i) = max(y(cuePoint:end));
        xlabel(['channel avg = ' num2str(channelMax(j,i))])
        ylabel(cue{j});
        set(gca,'FontSize',14)
        if j == 1
            title(['Electrode number: ' num2str(i)])
        end
        hold off
    end
    pause()
    close
end

disp('Complete')

% result: use the visual results for main selection and maxes for
% confirmation