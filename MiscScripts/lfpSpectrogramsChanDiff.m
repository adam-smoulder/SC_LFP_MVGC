% made to plot LFPs and their spectrograms to see if 300hz activity is
% present (artifact from filtering)
% load original dataset before running

% select intarg vs outtarg data and also what trial num to use (0 for
% random)
inTargVal = 1;
if exist('channelsToUse','var') == 0
    channelsToUse = [1 3];
end
if exist('inTargVal','var') == 0
    inTargVal = 1;
end
if exist('derivative','var') == 0
    derivative = 0;
end
if exist('demeanTrialAvg','var') == 0
    demeanTrialAvg = 0;
end
if exist('cueString','var') == 0
    cueString = 'sacclfp';
end

DS=30; % downsample factor of input data
maxPlotF = 30;

maxSgramPlotMag = 0; % use higher val for DS3
nFFT = 8092*8; % use higher val for DS3
pointsPerEval = 300/DS;
maxLFPPlotMag = 0;

fs = 30000/DS; % 10000 for DS3
%windowSize = fs/5;
windowSize = 300;
nOverlap = windowSize-pointsPerEval;
dur = 1;
if strcmp(cueString,'targlfp')
    startTime = 400;
elseif strcmp(cueString,'sacclfp')
    startTime = 600;
else
    startTime = 0;
end


t = 1:fs+1;

%%% organize data:
%%%
% inTarg = 0 -> out of target, inTarg = 1 -> in target, otherwise -> all
switch inTargVal
    case 0
        dataToUse = data([data.inTarg] == 0);
    case 1
        dataToUse = data([data.inTarg] == 1);
    otherwise
        dataToUse = data;
end

% sort channels from default to align with indexes
sortedData = zeros(16, length(t), length(dataToUse)); % dims: channel, time, trial
for trial = 1:length(dataToUse)
    disp(['Sorting trial ' num2str(trial)])
    trialData = dataToUse(trial);
    for i=1:length(trialData.channelOrder)
        eval(strcat('sortedData(trialData.channelOrder(i),:,trial) = trialData.', cueString, 'mat(i,:);'));
    end
end

if demeanTrialAvg
    sortedData = sortedData - repmat(mean(sortedData,3),1,1,size(sortedData,3));
end

if derivative
    temp = diff(sortedData,1,2);
    clear sortedData;
    sortedData = temp;
end


%%% Get Sgrams
%%% 
disp('Showing differences')

bipDiff = squeeze(mean(sortedData(channelsToUse(1),:,:),3))-squeeze(mean(sortedData(channelsToUse(2),:,:),3));
SDiff = abs(spectrogram(bipDiff, windowSize, nOverlap, nFFT, fs));
maxval = greatestMax(SDiff);
if maxLFPPlotMag == 0
    maxLFPPlotMag = greatestMax(abs(mean(sortedData(channelsToUse,:,:),3)));
end
if maxSgramPlotMag == 0
    maxSgramPlotMag = greatestMax(SDiff);
end



%%% Plotting
%%% 
plotTitle = 'No Time Difference Sgram';
if derivative
    plotTitle = '1st Difference Sgram';
end
if size(channelsToUse,1)==1
    yAxLab1 = ['Channels ' num2str(channelsToUse)];
else
    yAxLab1 = 'Channels a lot?';
end


sgramTime = 0:length(t)/(size(S,3)-1):length(t);
sgramFreqs = 0:(fs/2)/(size(S,2)-1):(fs/2);

figure(1)
subplot(1,2,1) %sgram
imagesc(sgramTime, sgramFreqs, SDiff)
colormap jet
set(gca, 'CLim', [0,maxSgramPlotMag]);
axis xy
axis([-inf inf 0 maxPlotF])
ylabel(yAxLab1);
colorbar
title(plotTitle)
hold on
plot(startTime*dur*ones(1,maxPlotF+1),0:maxPlotF,'r--','LineWidth',2)
hold off

subplot(1,2,2) %lfp
hold on
plot(squeeze(mean(sortedData(channelsToUse(1),:,:),3)));
plot(squeeze(mean(sortedData(channelsToUse(2),:,:),3)));
plot(bipDiff,'LineWidth',2);
axis([0 length(t)-1 -1.05*maxLFPPlotMag 1.05*maxLFPPlotMag])
axis xy
hold on
plot(startTime*dur*ones(1,2*ceil(maxLFPPlotMag)+1),-ceil(maxLFPPlotMag):ceil(maxLFPPlotMag),'r--','LineWidth',1)
hold off
legend('1st chan', '2nd chan','Difference','Location','NW')

