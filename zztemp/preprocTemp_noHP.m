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


% typically use [1 3; 5 7; 9 11] if none bad


%% PARAMETERS TO EDIT BEFORE RUNNING

fname = 'bl_sc_083017_mcell_spikelfp_cSC.mat';  % dataset file name
channelsToUse = [2 4; 6 8; 10 12];              % not "bad channels"
detrend = 2;                    % 0 = no detrend, 1 = pre-bip, 2 = post-bip
fs = 1000;                      % sampling rate of new data (Hz)
showAlignedSgrams = 0;          % show sgrams of unfiltered data
showDetrendSgrams = 0;          % show detrended sgrams
showBipSgrams = 0;              % show post-filtered sgrams
doNotch = 1;                    % 1 for do notch, 0 for no.
doHP = 0;                       % 1 for do hp, 0 for no.
harmonicNotch = 2;              % 0 = no, 1 = 5 harms, 2 = specific vals
inTargVal = 1;                  % intarg vs. outtarg vs. all data (used for visualization)

%% Load raw data mat

load(fname);                    % load data file
data = data([data.goodLFP]==1); % remove bad LFP trials
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
ChannelIsolationScript % rawData -> rawDataGoodChan
toc

%% Filtering
totalHPFilterOrder = 6;     % HP filter order (must be even # bc filtfilt)
highpFreq = 1;              % frequency for high pass filter
notchFreq = 60;             % base frequency to notch
notchOrder = 1;             % notch filter order (so 2*this+1 coeffs)

disp('Filtering data')
tic
FilteringScript % rawDataGoodChan -> lfp
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
    SgramFromAlignedData
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
for trial = 1:nTrials    
    % perform bipolar subtraction
    for j=1:length(channelsToUse) % for each pair of channels
        data(trial).targlfpbip(j,:) = data(trial).targlfpmat((2*j)-1,:) ...
            -data(trial).targlfpmat(2*j,:);
        data(trial).sacclfpbip(j,:) = data(trial).sacclfpmat((2*j)-1,:) ...
            -data(trial).sacclfpmat(2*j,:);
    end
    
%     % if 16 channeling...
%     data(trial).targlfpbip = data(trial).targlfpmat;
%     data(trial).sacclfpbip = data(trial).sacclfpmat;
end

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
else
    for i=1:length(data)
        % Store downsampled version of above as low SR lfp
        data(i).targlfpmat = data(i).targlfpbip;
        data(i).sacclfpmat = data(i).sacclfpbip;
    end
end

if showDetrendSgrams
    figureCount = figureCount+2;
    figure(figureCount)
    channelsForSgram = [1 2 3];
    SgramFromAlignedData
    hold on
    subplot(3,3,1)
    title(['Post - detrend, trial ' num2str(trialNum)])
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

%% SAVE data for use in GC
disp('Saving data')
tic

% remove fields we don't need that take up a lot of space
data = rmfield(data,[{'rawData'},{'targspikemat'},{'gospikemat'},...
    {'saccspikemat'},{'targlfpbip'},{'sacclfpbip'},{'golfpmat'},{'lfp'},...
    {'gazePosition'},{'targCode'},{'measCode'},{'goCode'},...
    {'stateTransitions'},{'srts'},{'delays'},{'spikeTimestamps'},...
    {'rawDataGoodChan'}]);
outputFileName = strcat('',fname(1:end-23),'_preproc_det_',...
    num2str(detrend),'_HP_',num2str(doHP),'_Notch_',num2str(harmonicNotch),...
    '_',num2str(fs),'hz.mat');

clearvars -except data detrend doHP doNotch DS figureCount fname ...
    fs highpFreq nnotch notchfreqs outputFileName rawFs totalHPFilterOrder ...
    channelsToUse

save(outputFileName,'-v7.3')

toc
disp('End Data Preparation')
