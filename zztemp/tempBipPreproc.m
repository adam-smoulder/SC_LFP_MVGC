% BAD channels by dataset (do NOT use these):
%
% bb 070915: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bb 071215: 11, 14 -> use [1 3; 6 8; 10 12]
% bb 080415: none
%
% bb 031315:  11 14 (use [1 3; 5 7; 9 12])
% 24bb 012317: none
% bb 121415: none
% bb 082917: none

% bl 071415: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 0723151: none (7, 11 is kinda bad) using [1 3; 5 8; 9 12]
% bl 0723152: none (7, 11 is kinda bad too) using [1 3; 5 8; 9 12]
%
% bl 031115: 2, 7, 11, 14
% bl 112515:  none
% bl 083017:  1

% typically use [1 3; 5 7; 9 11] if none bad


% channelsToUse = [1 3; 4 6; 8 10];              % not "bad channels"

%% PARAMETERS TO EDIT BEFORE RUNNING

supChans = [1 2 3 4];
intChans = [5 6 7 8];
deepChans = [9 10 11 12];
channelsToUse = 1:16;           % Fake; just the above variables are used
detrend = 2;                    % 0 = no detrend, 1 = pre-bip, 2 = post-bip
fs = 1000;                      % sampling rate of new data (Hz)
showAlignedSgrams = 0;          % show sgrams of unfiltered data
showDetrendSgrams = 0;          % show detrended sgrams
showBipSgrams = 0;              % show post-filtered sgrams
doNotch = 1;                    % 1 for do notch, 0 for no.
doHP = 1;                       % 1 for do hp, 0 for no.
harmonicNotch = 2;              % 0 = no, 1 = 5 harms, 2 = specific vals
inTargVal = 1;                  % intarg vs. outtarg vs. all data (used for visualization)

%% Load raw data mat

%load(fname);                    % load data file
data = data([data.goodLFP]==1); % remove bad LFP trials
%data = data([data.inTarg]==1);  % lets just use intarg
rawFs = 30000;                  % original sampling rate
DS = rawFs/fs;                  % downsample factor
figureCount = 0;                % counter for number of figures

if isfield(data,'lfp')
    data = rmfield(data,'lfp');
end

%% Isolate desired channels
disp(['Isolating desired channels: '...
    num2str(reshape(channelsToUse',[1,numel(channelsToUse)]))])

flatChanIndex = data(1).channelOrder;

for i=1:length(data)
       data(i).rawDataGoodChan = ...
           double(data(i).rawData(flatChanIndex,:));
end

tic

%% Filtering
totalHPFilterOrder = 6;     % HP filter order (must be even # bc filtfilt)
highpFreq = 1;              % frequency for high pass filter
notchFreq = 60;             % base frequency to notch
notchOrder = 1;             % notch filter order (so 2*this+1 coeffs)
rawFs = 30000;

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
    notchfreqs = [60, 161, 180, 300, 321.8, 483 14276, 14919.6]/(rawFs/2);
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
RealignScript
toc

if showAlignedSgrams
    figureCount = figureCount+1;
    figure(figureCount)
    SgramFromAlignedData16
end

%% Detrend (pre-BiP)

if detrend == 1
    disp('Detrending pre-Bipolar subtraction')
    tic
    PreBipDetrendScript
    toc
end


%% Bipolar subtraction
% subtract slightly separated channels within same general area to remove
% reference
% sacclfpmat -> sacclfpbip
disp('Performing Bipolar subtraction')
tic

nTrials = length(data);

PermutedBipolarSubtraction

if showBipSgrams
    figureCount = figureCount+1;
    figure(figureCount)
    SgramFromBipData
    hold on
    subplot(3,2,1)
    title(['Bipolar subtraction, trial ' num2str(trialNum)])
end

toc

%% Detrend (post-BiP)

if detrend == 2
    disp('Detrending post-Bipolar subtraction')
    tic
    PostBipDetrendScript
    toc
    if showDetrendSgrams
        figureCount = figureCount+2;
        figure(figureCount)
        %channelsForSgram = [1:16];
        SgramFromDetrendedData16
        hold on
        subplot(3,2,1)
        title(['Post - detrend, trial ' num2str(trialNum)])
    end
else
    for i=1:length(data)
        % Store downsampled version of above as low SR lfp
        data(i).targlfpmat = data(i).targlfpbip;
        data(i).sacclfpmat = data(i).sacclfpbip;
    end
end



%% Downsampling
if length(data(1).targlfpmat)-1 ~= fs  % only downsample if not done yet!
    disp(['Downsampling to ' num2str(fs) ' Hz'])
    tic
    for i=1:length(data)
        % Store downsampled version of above as low SR lfp
        tempTarg = resample(data(i).targlfpmat',1,DS)';
        tempSacc = resample(data(i).sacclfpmat',1,DS)';
        data(i).targlfpmat = tempTarg(1:3,:); % hackish to get relevant stuff
        data(i).sacclfpmat = tempSacc(1:3,:); 
    end
    toc
end

%% SAVE data for use in GC
disp('Saving data')
tic

% remove fields we don't need that take up a lot of space
data = rmfield(data,[{'rawData'},{'targspikemat'},{'gospikemat'},...
    {'saccspikemat'},{'targlfpbip'},{'sacclfpbip'},{'golfpmat'},{'lfp'},...
    {'gazePosition'},{'targCode'},{'measCode'},{'goCode'},...
    {'stateTransitions'},{'srts'},{'delays'},{'spikeTimestamps'},...
    {'rawDataGoodChan'}]);
outputFileName = strcat('BipPerm_',fname(1:end-23),'_filted_det_',...
    num2str(detrend),'_HP_',num2str(doHP),'_Notch_',num2str(harmonicNotch),...
    '_',num2str(fs),'hz.mat');

clearvars -except cueString data detrend doHP doNotch DS figureCount fname ...
    fs highpFreq nnotch notchfreqs outputFileName rawFs totalHPFilterOrder ...
    channelsToUse supChans intChans deepChans namesOfTheDay lNames theNameIWant

save(outputFileName,'-v7.3')

toc
disp('End Data Preparation')
