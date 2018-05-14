
% used to make sure channels for neural data show activity
% load preprocessed data beforehand
nChannels = 16;
nTime = fs+1;
nrt = 2; % number of random trials to use

% isolate intarg data
dataInTarg = data([data.inTarg] == 1);
nTrials = length(dataInTarg);

% extract data & get channel means
targVals = nan(nChannels, nTime, nTrials);
targVals = reshape([dataInTarg.targlfpmat],size(targVals));
targMean = squeeze(mean(targVals,3));
targStart = 400;
saccVals = reshape([dataInTarg.sacclfpmat],size(targVals));
saccMean = squeeze(mean(saccVals,3));
saccStart = 600;

% pick some random trials
trialsToUse = randi(nTrials,[1, nrt]);
targRandomTrials = targVals(:,:,trialsToUse);
saccRandomTrials = saccVals(:,:,trialsToUse);

figure()
% Plot trials and means
for i = 1:nChannels
    for j = 1:nrt
        subplot(nChannels, 2*nrt+2, (((2*nrt+2)*(i-1))+(2*j-1)))
        hold on
        %plot(squeeze(targRandomTrials(i,:,j)))
        %toPlot = conv(squeeze(targRandomTrials(i,:,j)),gausswin(15)); % smooths
        %plot(toPlot); %smooths
        %plot(targStart*ones(1,1000),-50:0.1:49.9,'r--','LineWidth',2)
        toPlot = squeeze(targRandomTrials(i,:,j));
        plot(-500:500,abs(fftshift(fft(toPlot))))
        axis([0 500 -inf inf])
        xlabel(['targ rando' num2str(j)])
        hold off
        
        subplot(nChannels, 2*nrt+2, (((2*nrt+2)*(i-1))+(2*j)))
        hold on
        toPlot = squeeze(saccRandomTrials(i,:,j));
        plot(-500:500,abs(fftshift(fft(toPlot))))
        axis([0 500 -inf inf])
        xlabel(['targ rando' num2str(j)])
        hold off
    end
    
    subplot(nChannels, 2*nrt+2, ((2*nrt+2)*(i-1))+2*nrt+1)
    hold on
    plot(squeeze(targMean(i,:)),'k.-','LineWidth',1.5)
    plot(targStart*ones(1,1000),-50:0.1:49.9,'r--','LineWidth',2)
    xlabel('targ avg')
    axis([0 nTime-1 min(targMean(i,:)) max(targMean(i,:))])
    hold off
    
    subplot(nChannels, 2*nrt+2, ((2*nrt+2)*(i-1))+2*nrt+2)
    hold on
    plot(squeeze(saccMean(i,:)),'b.-','LineWidth',1.5)
    plot(saccStart*ones(1,1000),-50:0.1:49.9,'r--','LineWidth',2)
    xlabel('sacc avg')
    axis([0 nTime-1 min(saccMean(i,:)) max(saccMean(i,:))])
    hold off
    
    subplot(nChannels, 2*nrt+2, ((2*nrt+2)*(i-1))+1)
    hold on
    ylabel(['Chan ' num2str(i)])
    hold off

end

