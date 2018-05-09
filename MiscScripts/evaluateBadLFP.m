% Evaluates where bad LFPs live


% bb 070915: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bb 071215: 11, 14 -> use [1 3; 6 8; 10 12]
% bb 080415: none
% bl 071415: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 0723151: none (7 is bad, 11 is kinda bad) -> use [1 3; 5 8; 9 12]
% bl 0723152: none

% typically use [1 3; 5 7; 9 11] if none bad


% channelsToUse = [1 3; 4 6; 8 10];              % not "bad channels"

%% 

fname = 'bl_sc0723152_mcell_spikelfp_cSC.mat';  % dataset file name
channelsToUse = [1 3; 5 8; 9 12];               % not "bad channels"
detrend = 0;                    % 0 = no detrend, 1 = pre-bip, 2 = post-bip
fs = 1000;                      % sampling rate of new data (Hz)
showAlignedSgrams = 1;          % show sgrams of unfiltered data
showDetrendSgrams = 1;          % show detrended sgrams
showBipSgrams = 1;              % show post-filtered sgrams
doNotch = 1;                    % 1 for do notch, 0 for no.
doHP = 1;                       % 1 for do hp, 0 for no.
harmonicNotch = 2;              % 0 = no, 1 = 5 harms, 2 = specific vals
inTargVal = 1;                  % intarg vs. outtarg vs. all data (used for visualization)

% Load raw data mat

load(fname);                    % load data file
data = data([data.goodLFP]==1); % remove bad LFP trials
rawFs = 30000;                  % original sampling rate
DS = rawFs/fs;                  % downsample factor
figureCount = 0;                % counter for number of figures

if isfield(data,'lfp')
    data = rmfield(data,'lfp');
end

% Isolate desired channels
disp(['Isolating desired channels: '...
    num2str(reshape(channelsToUse',[1,numel(channelsToUse)]))])
tic
ChannelIsolationScript % rawData -> rawDataGoodChan
toc

% Filtering
totalHPFilterOrder = 8;     % HP filter order (must be even # bc filtfilt)
highpFreq = 1;              % frequency for high pass filter
notchFreq = 60;             % base frequency to notch
notchOrder = 1;             % notch filter order (so 2*this+1 coeffs)

disp('Filtering data')
tic
FilteringScript % rawDataGoodChan -> lfp
toc

% Realign new lfp to targ, go, sacc - only these three lfp fields along with
% 'lfp' itself are updated by the end of this process

disp('Realigning data')

tic
RealignScript
toc




%%
trialNum = 1; % arbitrary
if ~exist('channelsForSgram')
    channelsForSgram = reshape(channelIndex',[1 6]); % select some channels to plot
end
cueString = 'targlfp';      % select cue

maxSgramPlotMag = 0;        % defaults to actual max
maxLfpPlotMag = 0;          % defaults to actual max
maxPlotF = 100;             % max frequency to display
nFFT = 8092*8;              % length of DFT
pointsPerEval = 300;        % step for center of shifting window
        
windowSize = rawFs/8;       % size of window
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
lChanPlot = length(channelsForSgram); % how many channels we're plotting
lFreq = nFFT/2+1;                   % size of frequency output for sgram
lTime = ceil(((rawFs+1)-windowSize)/pointsPerEval); % size of time output for sgram
trialLfp = zeros(lChanPlot,length(data(trialNum).targlfpmat(1,:)));
S = zeros(lChanPlot,lFreq,lTime);

% find the sgrams
for i = 1:lChanPlot
    evalString = strcat('trialLfp(i,:) = data(', num2str(trialNum), ').', cueString, 'mat(i,:);');
    eval(evalString);
    S(i,:,:) = abs(spectrogram(squeeze(trialLfp(i,:)), windowSize, nOverlap, nFFT, rawFs));
    disp(['Calculated sgram # ' num2str(i)])
end


% isolate intarg data
dataInTarg = data([data.inTarg] == 1);
nt = length(dataInTarg);

% extract data & get channel means
traceAvg = nan(lChanPlot, rawFs+1, nt);
targData = reshape([dataInTarg.targlfpmat],size(traceAvg));
meanTraces = mean(targData,3);


% setup for plotting sgrams
lfpTime = linspace(0, 1000, rawFs+1);
sgramTime = linspace(windowSize/DS, max(tSgram), lTime);
sgramFreqs = linspace(0, rawFs/2, lFreq);
maxLfpVal = greatestMax(trialLfp);
maxSgramVal = greatestMax(S);
channelLabels = reshape(channelsToUse',[1,6]);
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


for q = 1:nt
    trialLfp = squeeze(targData(:,:,q));
    figure(1)
    for i = 1:lChanPlot
        % plot trial
        subplot(lChanPlot,2,(2*i-1))
        plot(lfpTime,trialLfp(i,:),'LineWidth',0.5);
        axis([min(lfpTime) max(lfpTime) -1*maxLfpPlotMag maxLfpPlotMag])
        axis xy
        if i == 1
            title([cueString ', intarg = ' num2str(inTargVal)])
        end
        hold on
        ylabel(['Chan ' num2str(channelLabels(i))]);
        plot(startTime*ones(size([-1*maxLfpPlotMag:maxLfpPlotMag])),...
            [-1*maxLfpPlotMag:maxLfpPlotMag],'r--','LineWidth',0.5)
        if i == 1
            title(['Trial ' num2str(q) ' LFP'])
        end
        plot(startTime*ones(1,100),1:100,'r--','LineWidth',0.5)
        hold off
        
        
        % plot LFP for trial
        subplot(lChanPlot,2,2*i)
        plot(lfpTime(3:end),abs(diff(diff(trialLfp(i,:)))),'LineWidth',0.5);
        %axis([min(tSgram) max(tSgram) -1*maxLfpPlotMag maxLfpPlotMag])
        axis([-inf inf 0 50])
        axis xy
        if i == 1
            title('abs(2nd Derivative)')
        end
        hold on
        plot(startTime*ones(size([-1*maxLfpPlotMag:maxLfpPlotMag])),...
            [-1*maxLfpPlotMag:maxLfpPlotMag],'r--','LineWidth',0.5)
        plot(lfpTime(2:end),35*ones(length(lfpTime(2:end)),1))
        hold off
        
        %     % plot LFP mean trace
        %     subplot(lChanPlot,3,3*i)
        %     plot(lfpTime,meanTraces(i,:),'LineWidth',2);
        %     axis([min(tSgram) max(tSgram) -1*maxLfpPlotMag maxLfpPlotMag])
        %     axis xy
        %     if i == 1
        %         title([cueString ' mean, intarg = ' num2str(inTargVal)])
        %     end
        %     hold on
        %     plot(startTime*ones(size([-1*maxLfpPlotMag:maxLfpPlotMag])),...
        %         [-1*maxLfpPlotMag:maxLfpPlotMag],'r--','LineWidth',2)
        %     hold off
        
    end
    pause
end