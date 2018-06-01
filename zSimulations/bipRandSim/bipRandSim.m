%% Adapted from MVGC non-stationarity demo
%
% Demonstrates usage of the MVGC toolbox for non-stationary time series by
% "vertical" (windowed) regression on multi-trial data, for a minimal 2-variable
% VAR with linear trend and sinusoidally varying causal coefficient. Data is
% generated using the <var_to_tsdata_nonstat.html |var_to_tsdata_nonstat|>
% routine.
%
%% References
%
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014
% [ <matlab:open('mvgc_preprint.pdf') preprint> ].
%
% (C) Lionel Barnett and Anil K. Seth, 2012. See file license.txt in
% installation directory for licensing terms.
%
%% Parameters
tic
clear
useBip = 0;  %1 = bip sub, 2 = bip add
detrendIt = 2; % 1 = car before detrend, 2 = car after detrend
makeSomeNoise = 0; % 1 = WGN (Q), 2 = linear non-zero mean 2-step (R)
bipRand = 0; % do bipolar random subtraction
doICA = 0;
icaVarThresh = 0.1; % maximum variance for reference signal from ICA
noisePower = 1; % "SNR"
modelToUse = 4; % CHANGE MODEL HERE
shuffleCount = 5;
numruns = 10;


dnobs     = 0;       % initial observations to discard per trial - 0
regmode   = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS' ?
nlags     = 50;       % number of lags in VAR model
ntrials   = 100;    % number of trials - 1000
fs        = 500;%1000;     % sampling frequency (Hz) - 200
fres      = 4000;     % frequency resolution - 300?
nvars     = 3;%12;       % number of channels data that's being compared - 2
pointsPerEval = 20;  % evaluate GC at every ev-th sample - 5
modelOrder = 3;%17;       % model order (for real-world data should be estimated) - 2
dur       = 2; %0.1+modelOrder/fs;    % duration of trial (s) - 4.5
windowSize= dur*fs-modelOrder; % single window
seed      = 0;       % random seed (0 for unseeded) - 0

nobs  = dur*fs;                     % number of observations per trial
tnobs = nobs+dnobs;                 % total observations per trial for time series generation
k = 1:tnobs;                        % vector of time steps
t = (k-dnobs-1)/fs;                 % vector of times

%% Construct multiple non-stationary VAR time series

rng_seed(seed);

eval(strcat('model',num2str(modelToUse),'gen'))

trueX = X;

%% RANDOM BIP STUFF SETUP HERE

if bipRand == 1 % bip rand
    for datCount = 1:numruns
        clear X
        X = zeros(size(trueX)-[13,0,0]);
        refChans = 13:16;
        topChans = 1:4;
        midChans = 5:8;
        deepChans = 9:12;
    
        % for each trial, randbip it
        for i = 1:ntrials
            disp(['Bip replacement ' num2str(i) ' of ' num2str(ntrials)])
            
            % make reference for removal using random weights
            randoCoeffs = rand(length(refChans),1);
            coeffVec = randoCoeffs/sum(randoCoeffs);
            multMat = repmat(coeffVec,[1 nobs]);
            ref = squeeze(sum(multMat.*trueX(refChans,:,i)));
            
            % subtract the reference and reassign
            topTarg = squeeze(mean(trueX(topChans,:,i)))-ref;
            midTarg = squeeze(mean(trueX(midChans,:,i)))-ref;
            deepTarg = squeeze(mean(trueX(deepChans,:,i)))-ref;
            X(:,:,i) = [topTarg; midTarg; deepTarg];
        end
                
        
        % "Vertical" regression GC calculation
        nvars = size(X,1);
        
        wnobs = modelOrder+windowSize;                            % number of observations in "vertical slice"
        evalPoints    = wnobs : pointsPerEval : nobs;   % GC evaluation points
        enobs = length(evalPoints);                     % number of GC evaluations
        f1 = zeros(enobs,nvars,nvars,fres);               % (time, eq1, eq2, freq)
        nBins = fres+1;                                 % found mostly experimentally...
        stepsize = (fs/2)/nBins;                        % max frequency is fs/2 (Nyquist)
        freqs = 0:stepsize:(fs/2)-stepsize;
        
        % loop through evaluation points
        for e = 1:enobs
            j = evalPoints(e);
            fprintf('window %d of %d at time = %d',e,enobs,j);
            
            [AT,SIG] = tsdata_to_var(X(:,j-wnobs+1:j,:),modelOrder,regmode);
            if isbad(AT)
                fprintf(2,' *** skipping - VAR estimation failed\n');
                continue
            end
            
            [G,info] = var_to_autocov(AT,SIG);
            if info.error
                fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
                continue
            end
            if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
                fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
            end
            
            FF = autocov_to_pwcgc(G);
            if isbad(FF,false)
                fprintf(2,' *** skipping - GC calculation failed\n');
                continue
            end
            
            %Spectral GC against time
            a  = autocov_to_spwcgc(G,fres);
            dftLength1 = size(a,3);
            f1(e,:,:,1:dftLength1) = autocov_to_spwcgc(G,fres);
            fprintf('\n');
        end
        toc
        
        specGC = f1;
        timeGC = reshape(squeeze(sum(f1,4)),[size(f1,2),size(f1,3),size(f1,1)]);
        offset = modelOrder + windowSize;
        %specTime = t(offset:(length(t)-windowSize/2));
        specTime = t(offset:end);
        
        dispNullDists = 0;
        startTime = 0;
        inTargVal = 0;
        nullDistScript_trialShuffle;
        nullDistGC = squeeze(mean(specGC_perm))+2*squeeze(std(specGC_perm,0,1));
        
        specGC = max(0,real(specGC));
        
        
        if datCount == 1 % make the containers
            specGC_cont = zeros([numruns size(specGC)]);
            nullDistGC_cont = zeros([numruns size(nullDistGC)]);
        end
        specGC_cont(datCount,:,:,:,:) = specGC;
        nullDistGC_cont(datCount,:,:,:,:) = nullDistGC;
        
    end
    specGC = squeeze(mean(specGC_cont));
    nullDistGC = squeeze(mean(nullDistGC_cont));
    
else
    if bipRand == 2 % bip sub
        top = trueX(1,:,:)-trueX(3,:,:);
        mid = trueX(5,:,:)-trueX(7,:,:);
        deep = trueX(9,:,:)-trueX(11,:,:);
    else % just lumping if biprand isn't 1 or 2
        top = mean(trueX(1:4,:,:));
        mid = mean(trueX(5:8,:,:));
        deep = mean(trueX(9:12,:,:));
    end
    
    X = [top;mid;deep];
    
    % "Vertical" regression GC calculation
    nvars = size(X,1);
    
    wnobs = modelOrder+windowSize;                            % number of observations in "vertical slice"
    evalPoints    = wnobs : pointsPerEval : nobs;   % GC evaluation points
    enobs = length(evalPoints);                     % number of GC evaluations
    f1 = zeros(enobs,nvars,nvars,fres);               % (time, eq1, eq2, freq)
    nBins = fres+1;                                 % found mostly experimentally...
    stepsize = (fs/2)/nBins;                        % max frequency is fs/2 (Nyquist)
    freqs = 0:stepsize:(fs/2)-stepsize;
    
    % loop through evaluation points
    for e = 1:enobs
        j = evalPoints(e);
        fprintf('window %d of %d at time = %d',e,enobs,j);
        
        [AT,SIG] = tsdata_to_var(X(:,j-wnobs+1:j,:),modelOrder,regmode);
        if isbad(AT)
            fprintf(2,' *** skipping - VAR estimation failed\n');
            continue
        end
        
        [G,info] = var_to_autocov(AT,SIG);
        if info.error
            fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
            continue
        end
        if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
            fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
        end
        
        FF = autocov_to_pwcgc(G);
        if isbad(FF,false)
            fprintf(2,' *** skipping - GC calculation failed\n');
            continue
        end
        
        %Spectral GC against time
        a  = autocov_to_spwcgc(G,fres);
        dftLength1 = size(a,3);
        f1(e,:,:,1:dftLength1) = autocov_to_spwcgc(G,fres);
        fprintf('\n');
    end
    toc
    
    specGC = f1;
    timeGC = reshape(squeeze(sum(f1,4)),[size(f1,2),size(f1,3),size(f1,1)]);
    offset = modelOrder + windowSize;
    %specTime = t(offset:(length(t)-windowSize/2));
    specTime = t(offset:end);
    
    dispNullDists = 0;
    startTime = 0;
    inTargVal = 0;
    nullDistScript_trialShuffle;
    nullDistGC = squeeze(mean(specGC_perm))+2*squeeze(std(specGC_perm,0,1));
    
    specGC = max(0,real(specGC));
end

%% frequency plot 

if length(size(specGC)) > 3
    specGCShrunk = squeeze(sum(specGC,1));
else
    specGCShrunk = specGC;
end

if length(size(nullDistGC)) > 3
    nullDistShrunk = squeeze(sum(nullDistGC,1));
else
    nullDistShrunk = nullDistGC;
end

if size(specGCShrunk,3)~=length(freqs)
    freqs = freqs(1:size(specGCShrunk,3)); % hack, really don't feel like finding this rn...
end

numVar = size(specGC,2);
maxGC = greatestMax(specGCShrunk);

rows = 1:numVar;
cols = 1:numVar;

figure()
for i=1:length(rows)
    for j=1:length(cols)
        if rows(i)~=cols(j)
            subplot(numVar,numVar,(i-1)*numVar+j);
            plot(freqs,squeeze(specGCShrunk(rows(i),cols(j),:)))
            hold on
            plot(freqs,squeeze(nullDistShrunk(rows(i),cols(j),:)),'LineWidth',2)
            axis([0 250 0 maxGC])
        end
    end
end
subplot(numVar,numVar,1)
title(['row: ' num2str(min(rows)) '-' num2str(max(rows)) ' col: ' num2str(min(cols)) '-' num2str(max(cols))])
ylabel('GC Power')
xlabel('Frequency')

% for reference in figure
subplot(numVar,numVar,2)
legend('Actual','Null')
subplot(numVar, numVar, numVar^2)
xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
ylabel(['Model order = ' num2str(modelOrder)])


% establish the significance image and plot it
crunchedGC = squeeze(sum(specGCShrunk,3))/500;
crunchedNull = squeeze(sum(nullDistShrunk,3))/500;
sigImage = (crunchedGC > crunchedNull);

figure
subplot(1,3,1)
imagesc(crunchedGC)
title('GC')
set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
colorbar
subplot(1,3,2)
imagesc(crunchedNull)
title('Null Dist')
set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
colorbar
subplot(1,3,3)
imagesc(sigImage)
title('GCs above null dist')


%%
% <mvgc_demo_nonstationary.html back to top>
