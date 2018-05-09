trialNum = 15;              % arbitrary trial #
channelsToPlot = [1 2 3]; % select some channels to plot
cueString = 'targlfp';      % select cue

maxSgramPlotMag = 0;        % defaults to actual max
maxLfpPlotMag = 0;          % defaults to actual max
maxPlotF = 100;             % max frequency to display
nFFT = 8092*8;              % length of DFT
pointsPerEval = 300;        % step for center of shifting window
        
windowSize = rawFs/5;       % size of window
nOverlap = windowSize-pointsPerEval;
tSgram = linspace(0,1000,rawFs+1);  % since duration of trial is 1s

% organize/extract data:
% inTarg = 0 -> out of target, inTarg = 1 -> in target, otherwise -> all
switch inTargVal
    case 0
        dataToUse = data([data.inTarg] == 0);
    case 1
        dataToUse = data([data.inTarg] == 1);
    otherwise
        dataToUse = data;
end

% set up the sgram...
lChanPlot = length(channelsToPlot); % how many channels we're plotting
lFreq = nFFT/2+1;                   % size of frequency output for sgram
lTime = ceil(((rawFs+1)-windowSize)/pointsPerEval); % size of time output for sgram
trialLfp = zeros(lChanPlot,length(data(trialNum).targlfpbip(1,:)));
S = zeros(lChanPlot,lFreq,lTime);

% find the sgrams
for i = 1:lChanPlot
    evalString = strcat('trialLfp(i,:) = data(', num2str(trialNum), ').', cueString, 'bip(i,:);');
    eval(evalString);
    S(i,:,:) = abs(spectrogram(squeeze(trialLfp(i,:)), windowSize, nOverlap, nFFT, rawFs));
    disp(['Calculated sgram # ' num2str(i)])
end


% setup for plotting sgrams
lfpTime = linspace(0, 1000, rawFs+1);
sgramTime = linspace(0, max(tSgram), lTime);
sgramFreqs = linspace(0, rawFs/2, lFreq);
maxLfpVal = greatestMax(trialLfp);
maxSgramVal = greatestMax(S);
if maxLfpPlotMag == 0
    maxLfpPlotMag = maxLfpVal;
end 
if maxSgramPlotMag == 0
    maxSgramPlotMag = maxSgramVal;
end
if strcmp(cueString,'targlfp')
    startTime = 400;
elseif strcmp(cueString,'sacclfp')
    startTime = 600;
else
    startTime = 0;
end

% plotting
for i = 1:lChanPlot
    % plot sgram
    subplot(lChanPlot,2,(2*i-1))
    imagesc(sgramTime, sgramFreqs, squeeze(S(i,:,:)))
    ylabel(['Chan ' num2str(channelsToPlot(i))]);
    colormap jet
    set(gca, 'CLim', [0,maxSgramPlotMag]);
    axis xy
    axis([-inf inf 0 maxPlotF])
    if i == 1
        title(['Trial ' num2str(trialNum) ' LFP'])
    end
    colorbar
    hold on
    plot(startTime*ones(1,100),1:100,'r--','LineWidth',2)
    hold off
    
    % plot LFP
    subplot(lChanPlot,2,2*i)
    plot(lfpTime,trialLfp(i,:),'LineWidth',2);
    axis([min(tSgram) max(tSgram) -1*maxLfpPlotMag maxLfpPlotMag])
    axis xy
    if i == 1
        title([cueString ', intarg = ' num2str(inTargVal)])
    end
    hold on
    plot(startTime*ones(size([-1*maxLfpPlotMag:maxLfpPlotMag])),...
        [-1*maxLfpPlotMag:maxLfpPlotMag],'r--','LineWidth',2)
    hold off
end