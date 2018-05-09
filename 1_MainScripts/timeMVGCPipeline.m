% DOES TIME DOMAIN MVGC ONLY!!! this sums over all frequencies

% load appropriate preprocessed data before running!

% pref MO:
% bb 070915: 62
% bb 071215: 62
% bb 080415: 21 (62 is good)
% bb 031315: 62
% 24bb 012317: 
% bb 121415: 21 (62 is good)
% bb 082917: 21 (62 is good)

% bl 071415: 21 (62 is good)
% bl 0723151: 62 or 21
% bl 0723152: 21 (62 is good)
% bl 031115: 21 (62 is good)
% bl 112515: 25 (62 is good)
% bl 083017:


holdoutFraction = 0;  % for holdout/ensemble testing

%% GC Prep variables that may change between runs

cueString = 'sacclfp';  % cue to use
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

%dataToUse = dataToUse(randperm(length(dataToUse))); % scramble order
holdoutLength = ceil(holdoutFraction*length(dataToUse))+1;
useLength = length(dataToUse)-holdoutLength;
dataToUse = dataToUse(holdoutLength:end); % only take subset

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

%% Setup for MVGC
% Performs Multivariate Granger Causality using the MVGC toolbox
% if different model order desired, use following line:
% modelOrder = 60; 

% Parameters that should not change between runs
regmode   = 'OLS';          % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS')
nlags     = 100;             % number of lags in VAR model - 50 for performance sake
fres      = 8000*fs/1000;   % frequency resolution - use ~8000 for 1KHz
nvars     = 3;              % number of channels data that's being compared - 3 (superficial, mid, deep)
pointsPerEval = ceil(1*fs/1000);    % evaluate GC at every ev-th sample - 10 is pretty fair, lower is greater resolution
windowSize = 10*fs/1000;   % observation regression window size - 100ish
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
                  
timeGC = zeros(nvars, nvars, enobs);              % dims:  eq1, eq2, time
recover = 1;

%% "Vertical" regression GC calculation
disp('Beginning GC  calculation')
tic

% loop through evaluation points (where actual GC calculation is performed)
% note: evalPoints represents the leading edge of the window
for e = recover:enobs
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
    %var_info(info,true); % report results (and bail out on error)
    if info.error
        fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
        badCalcs(e) = 1;
        continue
    end    
    if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
        fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
        badCalcs(e) = 1;
    end
    
    % get GC for window
    F = autocov_to_pwcgc(G);
    assert(~isbad(F,false),'GC calculation failed');
    timeGC(:,:,e) = F;

    fprintf('\n');
end

toc
 
% represent windows with center timepoint of window
% offset = modelOrder + windowSize/2;
% specTime = t(offset:(length(t)-windowSize/2));
% %timeTime = t(offset:length(t)/enobs:(length(t)-windowSize/2));

% represent windows with leading edge
offset = modelOrder + windowSize;
timeTime = linspace(offset/fs,max(t),enobs);
%% plot results

numVar = size(timeGC,2);
maxGC = greatestMax(timeGC);

%TD Plot
figure(67)
for i=1:numVar
    for j=1:numVar
        if i~=j            
            subplot(numVar,numVar,(i-1)*numVar+j);
            plot(timeTime,squeeze(timeGC(i,j,:)), 'LineWidth', 3)
            ylabel('Granger Causality')
            axis xy
            axis([0, 1, 0, 1.2*maxGC])
            hold on
            plot(startTime*dur*ones(1,100),linspace(0,1.2*maxGC,100),'r--','LineWidth',2)
            hold off  
        end
    end
end

% for reference in figure
subplot(numVar, numVar, numVar^2)
title(['intarg = ' num2str(inTargVal)])
xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])

toc

%% Save relevant GC data
disp('Saving...')

clearvars -except AIC BIC cueString dnobs doHP doNotch DS dur fname freqs fres...
    fs inTargVal maxGC modelOrder numVar outputFileName pointsPerEval specGC ...
    specTime startTime t timeGC timeTime windowSize X channelsToUse supChans...
    intChans deepChans trialNums

biplabel = '';
if strcmp(outputFileName(1:3),'Bip')
    biplabel = 'Bip_';
end

outputFileName2 = strcat(biplabel,fname(1:end-23),...
    '_',cueString,...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');

save(outputFileName2,'-v7.3')

disp('MVGC Complete!')
