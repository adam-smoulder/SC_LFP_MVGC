% made to plot LFPs and their spectrograms to see if 300hz activity is
% present (artifact from filtering)
% load original dataset before running

% select intarg vs outtarg data and also what trial num to use (0 for
% random)
inTargVal = 1;
trialNum = 30;
channelsToPlot = 2;
cmax = 6000;
nFFT = 8092;
fs = 1000;
windowSize = 100;
pointsPerEval = 5;
nOverlap = windowSize-pointsPerEval;
maxPlotF = 500;




cueString = 'sacclfp';   % select cue
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
    trialData = dataToUse(trial);
    for i=1:length(trialData.channelOrder)
        eval(strcat('sortedData(trialData.channelOrder(i),:,trial) = trialData.', cueString, 'mat(i,:);'));
    end
end

if ~trialNum
    trialNum = randi(length(dataToUse));
end



%for any 1 trial
S = zeros(size(sortedData,1),nFFT/2+1,ceil((length(t)-windowSize)/pointsPerEval));
for i = 1:size(sortedData,1)
    S(i,:,:) = abs(spectrogram(squeeze(sortedData(i,:,trialNum)), windowSize, nOverlap, nFFT, fs));
end

sgramTime = 0:length(t)/(size(S,3)-1):length(t);
sgramFreqs = 0:(fs/2)/(size(S,2)-1):(fs/2);

if ~channelsToPlot
    channelsToPlot = 1;
end

maxval = greatestMax(S(channelsToPlot,:,:));

figure(2)
numChan2Plot = length(channelsToPlot);
for i = 1:numChan2Plot
    subplot(numChan2Plot,2,(2*i-1))
    imagesc(sgramTime, sgramFreqs, squeeze(S(channelsToPlot(i),:,:)))
    ylabel(['Chan ' num2str(channelsToPlot(i))]);
    colormap jet
    set(gca, 'CLim', [0,cmax]);
    axis xy
    axis([-inf inf 0 maxPlotF])
    if i == 1
        title(['Trial ' num2str(trialNum) ' LFP'])
    end
    
    subplot(numChan2Plot,2,2*i)
    plot(squeeze(sortedData(channelsToPlot(i),:,trialNum)));
    axis([0 length(t)-1 -150 150])
    axis xy
    if i == 1
        title([cueString ', intarg = ' num2str(inTargVal)])
    end
end


%% Same thing now, but after the bipolar subtraction
channelsToUse = [1 3; 5 7; 9 11];
demeanTrialAvg = 0;

if demeanTrialAvg
    demeanedData(:,:,:) = sortedData - repmat(mean(sortedData,3),1,1,size(sortedData,3));
else
    demeanedData = sortedData;
end

% for testing if demeaning by trial avg has any effect:
% demeanedData = sortedData;

% isolate data desired to be analyzed (subtracting is removing bias from
% voltage reference)
% for superficial, mid, and deep SC
X = zeros(3,1001, length(dataToUse));  % dims: depth (channel), time, trial
X(1,:,:) = squeeze(demeanedData(channelsToUse(1,1),:,:) - demeanedData(channelsToUse(1,2),:,:));
X(2,:,:) = squeeze(demeanedData(channelsToUse(2,1),:,:) - demeanedData(channelsToUse(2,2),:,:));
X(3,:,:) = squeeze(demeanedData(channelsToUse(3,1),:,:) - demeanedData(channelsToUse(3,2),:,:));

%for any 1 trial
S2 = zeros(3, nFFT/2+1, ceil((length(t)-windowSize)/pointsPerEval));
for i = 1:3
    S2(i,:,:) = abs(spectrogram(squeeze(X(i,:,trialNum)), windowSize, nOverlap, nFFT, fs));
end

maxval2 = greatestMax(S2);

figure(2)
for i = 1:3
    subplot(3,2,(2*i-1))
    imagesc(sgramTime, sgramFreqs, squeeze(S2(i,:,:)))
    ylabel(['Chan ' num2str(channelsToPlot(i))]);
    colormap jet
    set(gca, 'CLim', [0,500]);
    axis xy
    if i == 1
        title(['Trial ' num2str(trialNum) ' LFP'])
    end
    colorbar
    
    subplot(3,2,2*i)
    plot(squeeze(X(i,:,trialNum)));
    axis([0 1000 -150 150])
    axis xy
    if i == 1
        title([cueString ', intarg = ' num2str(inTargVal)])
    end
end
