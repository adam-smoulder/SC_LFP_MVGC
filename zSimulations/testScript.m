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
makeSomeNoise = 1; % 1 = WGN (Q), 2 = linear non-zero mean 2-step (R)
doICA = 0;
icaVarThresh = 0.1; % maximum variance for reference signal from ICA
noisePower = 1; % "SNR"
modelToUse = 6;
shuffleCount = 3;

dnobs     = 0;       % initial observations to discard per trial - 0
regmode   = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS' ?
nlags     = 50;       % number of lags in VAR model
ntrials   = 100;    % number of trials - 1000
dur       = 0.114;     % duration of trial (s) - 4.5
fs        = 1000;     % sampling frequency (Hz) - 200
fres      = 4000;     % frequency resolution - 300?
nvars     = 12;       % number of channels data that's being compared - 2
pointsPerEval = 20;  % evaluate GC at every ev-th sample - 5
windowSize= 100;     % observation regression window size - 75
%modelOrder= 2;       % model order (for real-world data should be estimated) - 2
seed      = 0;       % random seed (0 for unseeded) - 0

nobs  = dur*fs;                     % number of observations per trial
tnobs = nobs+dnobs;                 % total observations per trial for time series generation
k = 1:tnobs;                        % vector of time steps
t = (k-dnobs-1)/fs;                 % vector of times

%% Construct multiple non-stationary VAR time series

rng_seed(seed);

eval(strcat('model',num2str(modelToUse),'gen'))

%% do ICA if desired
if doICA
    %     stuff = reshape(X, nvars, ntrials*nobs); % placing tiles horizontally
    %     %%backup
    % %     stuff = nan(nvars,ntrials*nobs);
    % %     for i = 1:ntrials
    % %         stuff(:,((i-1)*nobs+1):(i*nobs)) = X(:,:,i);
    % %     end
    varDims = cell([1, ntrials]);
    removedRef = zeros(size(X));
    for i = 1:ntrials
        [ICs, AT, W] = fastica(squeeze(X(:,:,i))); % perform ICA
        [varA, minVarDim] = min(var(AT));   % find dimension with smallest variance
        if varA < icaVarThresh % if under thresh, remove suspected ref.
            removedRef(:,:,i) = AT(:,minVarDim)*ICs(minVarDim,:);
        end
        disp(['Trial ' num2str(i) ' of ' num2str(ntrials) ' complete'])
        varDims(i) = {var(AT)};
    end
    X = X-removedRef;
end

%%
predetX = X; % for reference later

if detrendIt == 1 %
    X = homeDetrend(X);
end

%bipolar subtraction or other
if (useBip == 1)||(useBip==2)
    tempX = X;
    clear X
    X = zeros(ceil(size(tempX)./[2 1 1]));
    for i = 1:size(X,1)
        if 2*i <= size(X,1) % if we have two channels to subtract
            if useBip ==1 % bip sub
                X(i,:,:) = tempX(2*i-1,:,:)-tempX(2*i,:,:);
            else % useBip == 2
                X(i,:,:) = tempX(2*i-1,:,:)+tempX(2*i,:,:);
            end
        else % odd number of vars and last channel
            X(i,:,:) = tempX(2*i-1,:,:);
        end
    end
end

if detrendIt == 2
    X = homeDetrend(X);
end

%% odd or superset comparisons
% oldX = X;
% chanA = sum(oldX(4,:,:),1);
% chanB = sum(oldX([1 2 3],:,:),1);
% 
% X = [chanA;chanB];

%% Estimation of model order
moEstX = X;

momax = 20; % max model order for estimation
icregmode = 'OLS';

ptic('\n*** tsdata_to_infocrit\n');
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(moEstX,momax,icregmode);
ptoc('*** tsdata_to_infocrit took ');

% Plot information criteria.

figure(1); clf;
order = 1:momax;
subplot(1,3,1)
hold on
plot(order,AIC,'LineWidth',2)
plot(order,BIC,'LineWidth',2)
title('Model order estimation');
subplot(1,3,2)
hold on
plot(order(2:end),diff(AIC),'LineWidth',2)
plot(order(2:end),diff(BIC),'LineWidth',2)
legend('AIC','BIC')
title('Derivative of model order ests.');
hold off
subplot(1,3,3)
hold on
plot(order(3:end),diff(AIC,2),'LineWidth',2)
plot(order(3:end),diff(BIC,2),'LineWidth',2)
legend('AIC','BIC')
title('2nd deriv');
hold off

fprintf('\nbest model order (AIC) = %d\n',moAIC);
fprintf('best model order (BIC) = %d\n',moBIC);

modelOrder = moAIC;


%% "Vertical" regression GC calculation
modelOrder = 14;
nvars = size(X,1);

wnobs = modelOrder+windowSize;                            % number of observations in "vertical slice"
evalPoints    = wnobs : pointsPerEval : nobs;   % GC evaluation points
enobs = length(evalPoints);                     % number of GC evaluations
f1 = zeros(enobs,nvars,nvars,fres);               % (time, eq1, eq2, freq)
nBins = fres+1;                                 % found mostly experimentally...
stepsize = (fs/2)/nBins;                        % max frequency is fs/2 (Nyquist)
freqs = 0:stepsize:(fs/2)-stepsize;

F21 = zeros(enobs,1); % will be our time GC values
F12 = zeros(enobs,1); % same^

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

%     F21(e) = FF(1,2); % estimated GC 2 -> 1 in time domain
%     F12(e) = FF(2,1); % estimated GC 1 -> 2 in time domain
    
    %Spectral GC against time
    a  = autocov_to_spwcgc(G,fres);
    dftLength1 = size(a,3);
    f1(e,:,:,1:dftLength1) = autocov_to_spwcgc(G,fres);
    
%     %Ensure time domain and spectral domain values align
%     fprintf('\nchecking that frequency-domain GC integrates to time-domain GC... \n');
%     Fint = smvgc_to_mvgc(a); % integrate spectral MVGCs
%     mad = maxabs(FF-Fint);
%     madthreshold = 1e-5;
%     if mad < madthreshold
%         fprintf('maximum absolute difference OK: = %.2e (< %.2e)\n',mad,madthreshold);
%     else
%         fprintf(2,'WARNING: high maximum absolute difference = %e.2 (> %.2e)\n',mad,madthreshold);
%     end
    
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


% plot GCs
% 
% numVar = size(specGC,2);
% maxGC = greatestMax(specGC(:,:,:,1:fres/2));
% 
% %SD Plot
% figure()
% for i=1:numVar
%     for j=1:numVar
%         if i~=j
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             imagesc(specTime,freqs,squeeze(specGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
%             ylabel('Frequency (Hz)')
%             axis xy
%             axis([-inf inf 0 100])
%             %axis([-inf, inf, 0, 25])
%             colormap jet
%             set(gca, 'CLim', [0,maxGC]);       
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar,numVar,1)
% title('ActualDist')
% subplot(numVar, numVar, numVar^2)
% set(gca, 'CLim', [0,maxGC]);
% c = colorbar;
% c.Label.String = 'GC';
% colorbar
% xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
% ylabel(['Model order = ' num2str(modelOrder)])
% 
% 
% nullDistMax = greatestMax(nullDistGC(:,:,:,1:fres/2));
% 
% %SD Plot
% figure()
% for i=1:numVar
%     for j=1:numVar
%         if i~=j
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             imagesc(specTime,freqs,squeeze(nullDistGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
%             ylabel('Frequency (Hz)')
%             axis xy
%             axis([-inf inf 0 100])
%             %axis([-inf, inf, 0, 25])
%             colormap jet
%             set(gca, 'CLim', [0,maxGC]);       
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar,numVar,1)
% title('NullDist')
% subplot(numVar, numVar, numVar^2)
% set(gca, 'CLim', [0,maxGC]);
% c = colorbar;
% c.Label.String = 'GC';
% colorbar
% xlabel(['Max nullDistGC = ' num2str(nullDistMax) '      points/eval = ' num2str(pointsPerEval)])
% ylabel(['Model order = ' num2str(modelOrder)])


%% frequency plot 

specGCShrunk = squeeze(sum(specGC,1));
if length(size(nullDistGC)) > 3
    nullDistShrunk = squeeze(sum(nullDistGC,1));
else
    nullDistShrunk = nullDistGC;
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
            axis([0 100 0 maxGC])
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
