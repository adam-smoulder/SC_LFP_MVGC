% using new bip replacement

%
% BAD channels by dataset (do NOT use these):
%
% bb 070915: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bb 071215: 11, 14 -> use [1 3; 6 8; 10 12]
% bb 080415: none - > [1 3; 5 7; 9 11] 
% bb 031315:  11 14 (use [1 3; 5 7; 9 12])
% bb 121415: none [1 3; 5 7; 9 11]
% bb 082917: none [1 3; 5 7; 9 11]
% 24bb 012317: none
% 24bb 081617: none

% bl 071415: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 0723151: none ->  using [1 3; 5 7; 9 11]
% bl 0723152: none -> using [1 3; 5 7; 9 11]
% bl 031115: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 112515:  none [1 3; 5 7; 9 11]
% bl 083017:  1 -> use [2 4; 6 8; 10 12]
% 24bl 071717: none
% 24bl 071718: none
% 24bl 071719: none


% channelsToUse = [1 3; 4 6; 8 10];              % not "bad channels"

%% PARAMETERS TO EDIT BEFORE RUNNING

%fname = 'bb_sc_080415_mcell_spikelfp_cSC.mat';  % dataset file name
channelsToUse = 1:16;           % not "bad channels"
detrend = 2;                    % 0 = no detrend, 1 = pre-bip, 2 = post-bip
fs = 1000;                      % sampling rate of new data (Hz)
showAlignedSgrams = 1;          % show sgrams of unfiltered data
showDetrendSgrams = 1;          % show detrended sgrams
showBipSgrams = 1;              % show post-filtered sgrams
doNotch = 1;                    % 1 for do notch, 0 for no.
doHP = 1;                       % 1 for do hp, 0 for no.
harmonicNotch = 2;              % 0 = no, 1 = 5 harms, 2 = specific vals
inTargVal = 1;                  % intarg vs. outtarg vs. all data (used for visualization)

% default channels for lumping
if ~exist('topChans','var')
    topChans = 1:4;
end
if ~exist('midChans','var')
    midChans = 5:8;
end
if ~exist('deepChans','var')
    deepChans = 9:12;
end
if ~exist('refChans','var')
    refChans = 13:16;
end

%% Load raw data mat

load(fname);                    % load data file
data = data([data.goodLFP]==1); % remove bad LFP trials
data = data([data.inTarg]==1);  % lets just use intarg
rawFs = 30000;                  % original sampling rate
DS = rawFs/fs;                  % downsample factor
figureCount = 0;                % counter for number of figures

if isfield(data,'lfp')
    data = rmfield(data,'lfp');
end

%% Isolate desired channels
disp(['Isolating desired channels: '...
    num2str(reshape(channelsToUse',[1,numel(channelsToUse)]))])
tic
%% Channel isolation script
% takes rawData -> rawDataGoodChan for the data variable
% only pulls out the selected channels for analysis
% -load data variable before running

flatChanIndex = data(1).channelOrder;

% only select desired channels
for i=1:length(data)
    data(i).rawDataGoodChan = ...
        double(data(i).rawData(flatChanIndex,:));
end
toc

%% Filtering
totalHPFilterOrder = 6;     % HP filter order (must be even # bc filtfilt)
highpFreq = 1;              % frequency for high pass filter
notchFreq = 60;             % base frequency to notch
notchOrder = 1;             % notch filter order (so 2*this+1 coeffs)

PERFORMLP = 0;              % perform lowpass, used for CSD

disp('Filtering data')
tic

% Setup high pass filter frequency and coeffs
highpNorm = highpFreq/(rawFs/2);
[bHigh,aHigh] = butter(floor(totalHPFilterOrder/2),highpNorm,'high');

lowpNorm = 250/(rawFs/2);
[bLow,aLow] = butter(floor(6/2),lowpNorm,'low');


% Setup notch filter frequencies and coeffs
if harmonicNotch == 1
    notchfreqs = [notchFreq:notchFreq:5*notchFreq 321.7]/(rawFs/2);
elseif harmonicNotch == 2 % visually inspected!
    notchfreqs = [60, 300, 321.8, 14276, 14919.6]/(rawFs/2);
else
    notchfreqs = [notchFreq 321.7]/(rawFs/2);
end
nnotch = length(notchfreqs);
anotch = zeros(nnotch,2*notchOrder+1); bnotch = zeros(nnotch,2*notchOrder+1);
for i=1:nnotch
    [bnotch(i,:),anotch(i,:)] = iirnotch(notchfreqs(i),notchfreqs(i)/30);
end

% Filter raw data for each trial
tic
for i=1:length(data)
    disp(['trial ' num2str(i)])
    rawData = double(data(i).rawDataGoodChan);
    
    notchedData = rawData;
    
    % Apply notch and high pass filters
    if doNotch
        for n=nnotch:-1:1
            notchedData = filtfilt(bnotch(n,:),anotch(n,:),notchedData')';
        end
    end
    
    if doHP
        data(i).lfp = filtfilt(bHigh,aHigh,notchedData')';
    else
        data(i).lfp = notchedData;
    end
    
    if PERFORMLP
        data(i).lfp = filtfilt(bLow,aLow,data(i).lfp')';
    end
    
end

disp(['notched @ ' num2str(notchfreqs*rawFs/2)])
toc

%% Realign Data
% Realign new lfp to targ, go, sacc - only these three lfp fields along with
% 'lfp' itself are updated by the end of this process

disp('Realigning data')

tic
% aligns lfp data about target onset and saccade onset
% takes lfp -> sacclfpmat and targlfpmat

for i=1:length(data)
    disp(['trial ' num2str(i)])
    trialData = double(data(i).lfp);
    
    % get data for cue times w.r.t. start of trial
    targCode = data(i).targCode;
    goCode = data(i).goCode;
    measCode = data(i).measCode;
    srt = data(i).srts;
    currstates = data(i).stateTransitions;
    
    targtime = double(currstates(2,currstates(1,:)==targCode))*30;
    gotime = double(currstates(2,currstates(1,:)==goCode))*30;
    sacctime = (double(currstates(2,currstates(1,:)==measCode))+srt)*30;
    
    % realign data for target onset
    targtimevec = (data(i).targtimevec(1)*30):(data(i).targtimevec(2)*30);
    data(i).targlfpmat = trialData(:,floor(targtime+targtimevec));
    
    % realign data for go cue
    if(~isempty(gotime))
        gotimevec = data(i).gotimevec(1)*30:data(i).gotimevec(2)*30;
        data(i).golfpmat = trialData(:,floor(gotime+gotimevec));
    end
    
    % realign data for saccade onset
    sacctimevec = data(i).sacctimevec(1)*30:data(i).sacctimevec(2)*30;
    data(i).sacclfpmat = trialData(:,floor(sacctime+sacctimevec));
end

toc

meanTargLFPPostRealign = mean(reshape([data.targlfpmat],[16, 30001, length(data)]),3);
meanSaccLFPPostRealign = mean(reshape([data.sacclfpmat],[16, 30001, length(data)]),3);

%% random bipolar reference creation and subtraction
nobs = size(data(1).targlfpmat,2);
nTrials = length(data);

for i = 1:nTrials
    disp(['Bip replacement ' num2str(i) ' of ' num2str(length(data))])
    
    % make reference for removal using random weights
    randoCoeffs = rand(length(refChans),1);
    coeffVec = randoCoeffs/sum(randoCoeffs);
    multMat = repmat(coeffVec,[1 nobs]);
    refTarg = sum(multMat.*data(i).targlfpmat(refChans,:));
    refSacc = sum(multMat.*data(i).sacclfpmat(refChans,:));
    
    % subtract the reference and reassign
    topTarg = mean(data(i).targlfpmat(topChans,:))-refTarg;
    midTarg = mean(data(i).targlfpmat(midChans,:))-refTarg;
    deepTarg = mean(data(i).targlfpmat(deepChans,:))-refTarg;
    data(i).targlfpmat = [topTarg; midTarg; deepTarg];

    topSacc = mean(data(i).sacclfpmat(topChans,:))-refSacc;
    midSacc = mean(data(i).sacclfpmat(midChans,:))-refSacc;
    deepSacc = mean(data(i).sacclfpmat(deepChans,:))-refSacc;
    data(i).sacclfpmat = [topSacc; midSacc; deepSacc];
end


%% Detrend

if detrend
    disp('Detrending post-Bipolar subtraction')
    tic
    % detrends data after bipolar subtraction
    
    % extract and detrend data
    tempTargData = zeros([size(data(1).targlfpmat) nTrials]);
    tempSaccData = zeros([size(data(1).sacclfpmat) nTrials]);
    tempTargData = homeDetrend(reshape([data.targlfpmat],size(tempTargData)));
    tempSaccData = homeDetrend(reshape([data.sacclfpmat],size(tempSaccData)));
    
    % reassign data
    for i=1:length(data)
        data(i).targlfpmat = squeeze(tempTargData(:,:,i));
        data(i).sacclfpmat = squeeze(tempSaccData(:,:,i));
    end
    
    disp('Finished detrending')
    toc
end



%% Downsampling
if length(data(1).targlfpmat)-1 ~= fs  % only downsample if not done yet!
    disp(['Downsampling to ' num2str(fs) ' Hz'])
    tic
    for i=1:length(data)
        % Store downsampled version of above as low SR lfp
        data(i).targlfpmat = resample(data(i).targlfpmat',1,DS)';
        data(i).sacclfpmat = resample(data(i).sacclfpmat',1,DS)';
    end
    toc
end

% these should be 0
meanTargLFPPostDS = mean(reshape([data.targlfpmat],[3, 1001, length(data)]),3);
meanSaccLFPPostDS = mean(reshape([data.sacclfpmat],[3, 1001, length(data)]),3);

%% 
% 
% SAVE data for use in GC
disp('Saving data')
tic

% remove fields we don't need that take up a lot of space
data = rmfield(data,[{'rawData'},{'targspikemat'},{'gospikemat'},...
    {'saccspikemat'},{'golfpmat'},{'lfp'},...
    {'gazePosition'},{'targCode'},{'measCode'},{'goCode'},...
    {'stateTransitions'},{'srts'},{'delays'},{'spikeTimestamps'},...
    {'rawDataGoodChan'}]);

theTime = clock;
endPart = strcat(num2str(theTime(2)),num2str(theTime(3)),num2str(theTime(4)),num2str(theTime(5)));

outputFileName = strcat('preprocRandBip_',fname(1:end-23),'_',endPart);

clearvars -except data detrend doHP doNotch DS figureCount fname ...
    fs highpFreq nnotch notchfreqs outputFileName rawFs totalHPFilterOrder ...
    meanTargLFPPostRealign meanSaccLFPPostRealign endPart numruns datCount ...
    targSpecGC_holder saccSpecGC_holder targNullGC_holder saccNullGC_holder


save(outputFileName,'-v7.3')

toc
disp('End Data Preparation')

