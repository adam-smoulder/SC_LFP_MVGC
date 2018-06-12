% load data saved from MVGC run before using this!

% shows beta band power for event-related activity
% [500 900] for target onset (target appears at 400)
% [400 800] for saccade onset (saccade cue is at 600)

% decide time window
if strcmp(cueString,'targlfp')
    window = [300 900]/fs;
else
    window = [400 800]/fs;
end

specTime = linspace(specTime(1),specTime(end),size(specGC,1));

% integrate in time across time window
timesToUse = (specTime > window(1)) & (specTime < window(2)); % timeTime maps ms to the STFT windows
newSpecTimes = specTime(timesToUse);
convFact = (size(specGC,4)-1)/(fs/2); % max frequency = nyquist
betaBand = (12*convFact):(30*convFact); % 12-30 Hz
bbp = squeeze(sum(specGC(timesToUse,:,:,betaBand),4)); % betaBandPower: dims: eq x eq x freq

% plot
maxMagVal = greatestMax(bbp);
maxFreqToPlot = 100;

figure(123)

for i = 1:numVar
    for j = 1:numVar
        if j~=i
            subplot(numVar,numVar,(i-1)*numVar+j)
            plot(newSpecTimes,squeeze(bbp(:,i,j)),'LineWidth',3)
            axis([-inf, inf, 0, 1.2*maxMagVal])
            hold on
            plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
            hold off   
            set(gca,'fontsize',14)
        end
    end
end

subplot(numVar, numVar,2)
title(['Beta Band, ' cueString])
set(gca,'fontsize',14)
subplot(numVar,numVar,4)
ylabel('Granger Causality Magnitude')
set(gca,'fontsize',14)
subplot(numVar,numVar,8)
xlabel('Time (ms)')
set(gca,'fontsize',14)

