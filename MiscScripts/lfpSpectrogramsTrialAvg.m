% made to plot LFPs and their spectrograms to see if 300hz activity is
% present (artifact from filtering)
% load original dataset before running!

% select intarg vs outtarg data and also what trial num to use (0 for
% random)
if exist('inTargVal','var') == 0
    inTargVal = 1;
end
if exist('channelsToPlot','var') == 0
    channelsToPlot = [1 3]; % select some channels to plot
end
if exist('cueString','var') == 0
    cueString = 'sacclfp';   % select cue
end
if exist('derivative','var') == 0
    derivative = 0; % take derivative?
end
DS=30; % downsample factor of input data
demeanTrialAvg = 0;

maxPlotF = 100;

cmax = 0; % use higher val for DS3
nFFT = 8092*4; % use higher val for DS3
pointsPerEval = 10*30/DS;
maxPlotMag = 0;

fs = 30000/DS; % 10000 for DS3
%windowSize = fs/5;
windowSize = fs*3/10; % idk
nOverlap = windowSize-pointsPerEval;
dur = 1;
if strcmp(cueString,'targlfp')
    startTime = 400*30/DS;
elseif strcmp(cueString,'sacclfp')
    startTime = 600*30/DS;
else
    startTime = 0;
end


t = 1:fs+1;
% organize data:
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

if ~channelsToPlot
    channelsToPlot = 1;
end

%for nonbipolar stuff
disp(['Making trial sgrams'])
% sometimes have to delete the +1...
S = zeros(length(channelsToPlot),nFFT/2+1,ceil((length(t)-windowSize)/pointsPerEval));
for i = 1:length(channelsToPlot)
    S(i,:,:) = abs(spectrogram(squeeze(mean(sortedData(channelsToPlot(i),:,:),3)), windowSize, nOverlap, nFFT, fs));
end

sgramTime = 0:length(t)/(size(S,3)-1):length(t);
sgramFreqs = 0:(fs/2)/(size(S,2)-1):(fs/2);


maxval = greatestMax(S);
if maxPlotMag == 0
    maxPlotMag = greatestMax(abs(sortedData(channelsToPlot,:,:)));
end

figure(1)
numChan2Plot = length(channelsToPlot);
if cmax == 0
    cmax = maxval;
end

showSub=0;
if length(channelsToPlot)==2
    showSub = 1;
end

disp(['Plotting channel sgrams'])
for i = 1:numChan2Plot
    subplot(numChan2Plot+showSub,2,(2*i-1))
    imagesc(sgramTime, sgramFreqs, squeeze(S(i,:,:)))
    ylabel(['Chan ' num2str(channelsToPlot(i))]);
    colormap jet
    set(gca, 'CLim', [0,cmax]);
    axis xy
    axis([min(sgramTime) max(sgramTime) 0 maxPlotF])
    if i == 1
        title(['Trial Average LFP'])
    end
    colorbar
    hold on
    plot(startTime*dur*ones(1,maxPlotF+1),0:maxPlotF,'r--','LineWidth',2)
    hold off
    
    subplot(numChan2Plot+showSub,2,2*i)
    plot(squeeze(mean(sortedData(channelsToPlot(i),:,:),3)),'LineWidth',2);
    axis([0 length(t)-1 -1*maxPlotMag maxPlotMag])
    axis xy
    if i == 1
        title([cueString ', intarg = ' num2str(inTargVal)])
    end
    hold on
    plot(startTime*dur*ones(1,2*ceil(maxPlotMag)+1),-ceil(maxPlotMag):ceil(maxPlotMag),'r--','LineWidth',1) 
    hold off
end
if showSub
    disp(['Showing differences'])

    theDiff = squeeze(mean(sortedData(channelsToPlot(1),:,:),3))-squeeze(mean(sortedData(channelsToPlot(2),:,:),3));
    SDiff = abs(spectrogram(theDiff, windowSize, nOverlap, nFFT, fs));
    
    sgramTime = 0:length(t)/(size(S,3)-1):length(t);
    sgramFreqs = 0:(fs/2)/(size(S,2)-1):(fs/2);

    subplot(3,2,5)
    imagesc(sgramTime, sgramFreqs, SDiff)
    ylabel('Difference');
    colormap jet
    set(gca, 'CLim', [0,cmax]);
    axis xy
    axis([min(sgramTime) max(sgramTime) 0 maxPlotF])
    colorbar
    title('Difference')
    hold on
    plot(startTime*dur*ones(1,maxPlotF+1),0:maxPlotF,'r--','LineWidth',2)
    hold off
    
    subplot(3,2,6)
    hold on
    plot(squeeze(mean(sortedData(channelsToPlot(1),:,:),3)));
    plot(squeeze(mean(sortedData(channelsToPlot(2),:,:),3)));
    plot(theDiff,'LineWidth',2);
    axis([0 length(t)-1 -1*maxPlotMag maxPlotMag])
    axis xy
    title('Difference')
    hold on
    plot(startTime*dur*ones(1,2*ceil(maxPlotMag)+1),-ceil(maxPlotMag):ceil(maxPlotMag),'r--','LineWidth',1) 
    hold off
    legend('1st chan', '2nd chan','Difference','Location','NW')

end
