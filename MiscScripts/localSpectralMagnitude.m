% load data saved from MVGC run before using this!

% the purpose of this script is to show the magnitude of the spectral GC
% within the window of activity of a given trial, which corresponds to:
% [500 900] for target onset (target appears at 400)
% [400 800] for saccade onset (saccade cue is at 600)

% decide time window
if strcmp(cueString,'targlfp')
    window = [400 800]/fs;
else
    window = [400 800]/fs;
end

specTime = linspace(specTime(1),specTime(end),size(specGC,1));

% integrate in time across time window
timesToUse = (specTime > window(1)) & (specTime < window(2)); % timeTime maps ms to the STFT windows
newSpecTimes = specTime(timesToUse);
localSpecGCMag = squeeze(sum(specGC(timesToUse,:,:,:))); % dims: eq x eq x freq

% plot
maxMagVal = greatestMax(localSpecGCMag);
maxFreqToPlot = 100;

figure(123)

for i = 1:numVar
    for j = 1:numVar
        if j~=i
            subplot(numVar,numVar,(i-1)*numVar+j)
            plot(freqs,squeeze(localSpecGCMag(i,j,:)),'LineWidth',3)
            axis([0, maxFreqToPlot, 0, 1.2*maxMagVal])
            
            set(gca,'fontsize',14)
        end
    end
end

subplot(numVar, numVar,2)
title('Saccade Onset -> 400ms after')
set(gca,'fontsize',14)
subplot(numVar,numVar,4)
ylabel('Granger Causality Magnitude')
set(gca,'fontsize',14)
subplot(numVar,numVar,8)
xlabel('Frequency (Hz)')
set(gca,'fontsize',14)

