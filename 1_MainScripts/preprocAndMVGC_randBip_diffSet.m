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

fname = 'bb_sc_080415_mcell_spikelfp_cSC.mat';  % dataset file name
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

%% bipolar replacement
nobs = size(data(1).targlfpmat,2);
nTrials = length(data);

for i = 1:nTrials
    disp(['Bip replacement ' num2str(i) ' of ' num2str(length(data))])
    randoCoeffs = rand(4,1);
    coeffVec = randoCoeffs/sum(randoCoeffs);
    multMat = repmat(coeffVec,[1 nobs]);
    
    % make reference for removal
    refTarg = sum(multMat.*data(i).targlfpmat(refChans,:));
    refSacc = sum(multMat.*data(i).sacclfpmat(refChans,:));
    
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
    meanTargLFPPostRealign meanSaccLFPPostRealign endPart

save(outputFileName,'-v7.3')

toc
disp('End Data Preparation')
%%






%%%
% 
% RUNNING MVGC
% 
%%%










%%%
%% GC Prep variables that may change between runs

cueString = 'targlfp';  % cue to use
inTargVal = 1;          % intarg vs outtarg vs both
analysisType = 'cond';  % 'cond' = conditional, 'pw' = pairwise

holdoutFraction = 0;  % for holdout/ensemble testingl; 0 for using all data
nullDistChoice = 'trial';   % trial for trial scramble, time for time scramble


% Setup for GC
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
momax     = 70;     % maximum model order for model order estimation

% Parameters that should not change between runs
dnobs     = 0;          % initial observations to discard per trial - 0; we can ignore as needed
dur       = 1;          % duration of trial (s) - 1s

nobs  = dur*fs+1;       % number of observations in a trial
tnobs = nobs+dnobs;     % total observations per trial for time series generation
k = 1:tnobs;            % vector of time steps

monkey_data = fname(1:12);          % used for file naming and such
axisLimit = 0;                      % if 0, axisLimit = maxGC (only used for subtractorExtractor)
downsampleFactor = 1000/fs;         % HACK fix later

%% Isolate desired data

% inTarg = 0 -> out of target, inTarg = 1 -> in target, otherwise -> all
switch inTargVal
    case 0
        dataToUse = data([data.inTarg] == 0 & [data.goodLFP] == 1);
    case 1
        dataToUse = data([data.inTarg] == 1  & [data.goodLFP] == 1);
    otherwise
        dataToUse = data;
end

% remove held out data if desired
if holdoutFraction
    holdoutLength = ceil(holdoutFraction*length(dataToUse))+1;
    useLength = length(dataToUse)-holdoutLength;
    dataToUse = dataToUse(holdoutLength:end); % only take subset
end

% specifically isolate the cue we want to analyze
[nChannels, nTime] = size(dataToUse(1).sacclfpmat);
nTrials = length(dataToUse);
X = nan(nChannels, nTime, nTrials);
if strcmp(cueString,'targlfp')
    X = reshape([dataToUse.targlfpmat],size(X));
elseif strcmp(cueString,'sacclfp')
    X = reshape([dataToUse.sacclfpmat],size(X));
else
    disp('Something went wrong...Check your cueString')
end

trialNums = [data.trialNum];

% remove the rest of the stuff from memory, we shouldn't need it!
clear data;
clear dataToUse;


%% Model order estimation
X4mo = X;
% target
% 200-900 is normal
% 300-800 is like pre
% 200-800 is like pre
% 300-900 is like pre

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X4mo,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria. 
figureCount = figureCount+1;
figure(figureCount); clf;
order = 1:momax;
subplot(1,2,1)
hold on
plot(order,AIC,'LineWidth',2)
plot(order,BIC,'LineWidth',2)
title(['Model order estimation']);
legend('AIC','BIC')
subplot(1,2,2)
hold on
plot(order(2:end),diff(AIC),'LineWidth',2)
plot(order(2:end),diff(BIC),'LineWidth',2)
legend('AIC','BIC')
title('1st Difference of model order est.');
hold off

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

modelOrder = moAIC;
clear X4mo

%% Setup for MVGC
% Performs Multivariate Granger Causality using the MVGC toolbox
% if different model order desired, use following line:
% modelOrder = 62; 

% Parameters that should not change between runs
regmode   = 'OLS';          % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS')
nlags     = 100;             % number of lags in VAR model - 50 for performance sake
fres      = 8000*fs/1000;   % frequency resolution - use ~8000 for 1KHz
nvars     = 3;              % number of channels data that's being compared - 3 (superficial, mid, deep)
pointsPerEval = ceil(10*fs/1000);    % evaluate GC at every ev-th sample - 10 is pretty fair, lower is greater resolution
windowSize = 100*fs/1000;   % observation regression window size - 100ish
tstat     = 'F';            % 'F' for Granger's F-Test, 'chi2' for chi-squared
alpha     = 0.05;           % significance level
mhtc      = 'FDR';          % multiple hypothesis test correction (see routine 'significance')

if strcmp(cueString,'targlfp')
    startTime = 0.4;
elseif strcmp(cueString,'sacclfp')
    startTime = 0.6;
else
    startTime = 0;
end

t = (k-dnobs-1)/fs;                             % vector of times
wnobs = modelOrder+windowSize;                  % number of observations in "vertical slice"
evalPoints = wnobs:pointsPerEval:nobs;          % GC evaluation points
enobs = length(evalPoints);                     % number of GC evaluations
nBins = fres+1;                                 % found mostly experimentally...
stepsize = (fs/2)/nBins;                        % max frequency is fs/2 (Nyquist)
freqs = 0:stepsize:(fs/2)-stepsize;             % frequency values
badCalcs = zeros(1,enobs);                      % record where calculations fail
                  
specGC = zeros(enobs,nvars,nvars,nBins);        % dims:  time, eq1, eq2, freq
timeGC = zeros(nvars, nvars, enobs);            % dims:  eq1, eq2, time

%% "Vertical" regression GC calculation
disp('Beginning GC  calculation')
disp(['Time MVGC calculated up to ' num2str(fs/10) 'Hz'])
tic

% loop through evaluation points (where actual GC calculation is performed)
% note: evalPoints represents the leading edge of the window
for e = 1:enobs
    j = evalPoints(e);
    fprintf('window %d of %d at time = %d',e,enobs,j);
    
    % convert time series data in window to VAR
    [A,SIG] = tsdata_to_var(X(:,j-wnobs+1:j,:),modelOrder,regmode);
    if isbad(A)
        fprintf(2,' *** skipping - VAR estimation failed\n');
        badCalcs(e) = 1;
        continue
    end
    
    % calculate autocovariance from VAR
    [G,info] = var_to_autocov(A,SIG);
    if info.error
        fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
        badCalcs(e) = 1;
        continue
    end
    if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
        fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
        badCalcs(e) = 1;
    end
    
    
    % Calculate GC for window
    a  = autocov_to_spwcgc(G,fres); % rate limiting step
    specGC(e,:,:,1:size(a,3)) = a;
    timeGC(:,:,e) = squeeze(sum(specGC(e,:,:,1:floor(nBins/5)),4))./(nBins/5); % up to (nyq/5)
    if isbad(a,false)
        fprintf(2,' *** skipping - GC calculation failed\n');
        badCalcs(e) = 1;
        continue
    end

    fprintf('\n');
end

toc
 
% represent windows with center timepoint of window
% offset = modelOrder + windowSize/2;
% specTime = t(offset:(length(t)-windowSize/2));
% %timeTime = t(offset:length(t)/enobs:(length(t)-windowSize/2));

% represent windows with leading edge
offset = modelOrder + windowSize;
specTime = t(offset:end);

origSpecGC = specGC;
origTimeGC = timeGC;

% negative and imaginary values are uninterpretable, so:
% find real part, then make 0 for all neg. vals
specGC = max(0,real(specGC)); 
timeGC = max(0,real(timeGC));

%% Calculate null distribution and subtract it

shuffleCount = 10;
dispNullDists = 0;

nullDistChoice = ternaryOp(exist('nullDistChoice','var'),nullDistChoice,'trial'); % default it

% find null distribution (or set to 0 if not desired)
if strcmp(nullDistChoice,'trial')
    nullDistScript_trialShuffle % outputs specGC_perm and timeGC_perm
elseif strcmp(nullDistChoice,'time')
    nullDistScript_timeScramble % outputs specGC_perm and timeGC_perm
else 
    specGC_perm = zeros([shuffleCount size(specGC)]);
end

% update spec and time GC based on null dist
nullDistSpecGC = squeeze(mean(specGC_perm,1));
nullDistTimeGC = squeeze(mean(timeGC_perm,1));

%% plot results

numVar = size(specGC,2);
maxGC = greatestMax(specGC(:,:,:,1:fres/2));

% spectrogram
figure()
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(specGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
            ylabel('Frequency (Hz)')
            axis xy
            axis([0 1 0 50])
            %axis([-inf, inf, 0, 25])
            colormap jet
            set(gca, 'CLim', [0,maxGC]);
            hold on
            plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
            hold off            
        end
    end
end

% for reference in figure
subplot(numVar, numVar, numVar^2)
set(gca, 'CLim', [0,maxGC]);
c = colorbar;
c.Label.String = 'GC';
colorbar
title(['intarg = ' num2str(inTargVal)])
xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
ylabel(['Model order = ' num2str(modelOrder)])

subplot(numVar, numVar, 2)
title([fname(1:2) fname(7:12)])



maxNullGC = greatestMax(nullDistSpecGC(:,:,:,1:fres/2));

% spectrogram for null dist
figure()
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(nullDistSpecGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
            ylabel('Frequency (Hz)')
            axis xy
            axis([0 1 0 50])
            %axis([-inf, inf, 0, 25])
            colormap jet
            set(gca, 'CLim', [0,maxGC]);
            hold on
            plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
            hold off            
        end
    end
end

% for reference in figure
subplot(numVar, numVar, numVar^2)
set(gca, 'CLim', [0,maxGC]);
c = colorbar;
c.Label.String = 'GC';
colorbar
title(['intarg = ' num2str(inTargVal)])
xlabel(['Max GC = ' num2str(maxNullGC) '      points/eval = ' num2str(pointsPerEval)])
ylabel(['Model order = ' num2str(modelOrder)])

subplot(numVar, numVar, 2)
title([fname(1:2) fname(7:12)])

% %TD Plot
% maxTimeGC = greatestMax(timeGC);
% figure(67)
% for i=1:numVar
%     for j=1:numVar
%         if i~=j            
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             plot(timeTime,squeeze(timeGC(i,j,:)), 'LineWidth', 3)
%             ylabel('Granger Causality')
%             axis xy
%             axis([-inf, inf, 0, 1.2*maxTimeGC])
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar, numVar, numVar^2)
% title(['intarg = ' num2str(inTargVal)])
% xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])

totalTime = toc

%% Save relevant GC data
disp('Saving...')

clearvars -except AIC BIC cueString dnobs doHP doNotch DS dur fname freqs fres...
    fs inTargVal maxGC modelOrder numVar outputFileName pointsPerEval specGC ...
    specTime startTime t timeGC timeTime windowSize X channelsToUse supChans...
    intChans deepChans trialNums origSpecGC origTimeGC specGC_perm timeGC_perm ...
    nullDistChoice endPart totalTime

pretag = ['randBip' endPart '_'];

outputFileName2 = strcat(pretag,fname(1:end-23),...
    '_',cueString,...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');

save(outputFileName2,'-v7.3')

disp('MVGC Complete!')

