% load appropriate preprocessed data before running! Change range of times
% to use of Model order testing below

% GC Prep variables that may change between runs

cueString = 'targlfp';  % cue to use
inTargVal = 1;          % intarg vs outtarg vs both
analysisType = 'cond';  % 'cond' = conditional, 'pw' = pairwise

% Setup for GC
icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
momax     = 70;     % maximum model order for model order estimation

% Parameters that should not change between runs
dnobs     = 0;          % initial observations to discard per trial - 0; we can ignore as needed
dur       = 1;          % duration of trial (s) - 1s

nobs  = dur*fs+1;       % number of observations in a trial
tnobs = nobs+dnobs;     % total observations per trial for time series generation
k = 1:tnobs;            % vector of time steps

monkey_data = fname(1:12);          % 
axisLimit = 0;                      % if 0, axisLimit = maxGC (only used for subtractorExtractor)
downsampleFactor = 1000/fs;         % takes 1khz -> 250hz % HACK fix later


% Isolate desired data

% inTarg = 0 -> out of target, inTarg = 1 -> in target, otherwise -> all
switch inTargVal
    case 0
        dataToUse = data([data.inTarg] == 0 & [data.goodLFP] == 1);
    case 1
        dataToUse = data([data.inTarg] == 1  & [data.goodLFP] == 1);
    otherwise
        dataToUse = data;
end

dataToUse = dataToUse(randperm(length(dataToUse))); % scramble order


% specifically isolate the cue we want to analyze
[nChannels, nTime] = size(dataToUse(1).sacclfpmat);
nTrials = length(dataToUse);  % updated now for in vs. outtarg
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
range = [450 650];
%X4mo = X(:,range(1):range(2),:);
X4mo = X;

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X4mo,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria. 
%figureCount = figureCount+1;
figure(); clf;
order = 1:momax;
hold on
plot(order,AIC,'LineWidth',2)
plot(order,BIC,'LineWidth',2)
title(['Model order estimation, ' num2str(range)]);
legend('AIC','BIC')
hold off

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

modelOrder = moAIC;