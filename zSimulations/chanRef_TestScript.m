%% CHAN REF TESTING load the thing here:
load('X_forChanRefTests_predet');
% Parameters that should not change between runs
regmode   = 'OLS';          % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS')
nlags     = 100;             % number of lags in VAR model - 50 for performance sake
fres      = 8000*fs/1000;   % frequency resolution - use ~8000 for 1KHz
nvars     = 12;              % number of channels data that's being compared - 3 (superficial, mid, deep)
pointsPerEval = ceil(10*fs/1000);    % evaluate GC at every ev-th sample - 10 is pretty fair, lower is greater resolution
windowSize = 100*fs/1000;   % observation regression window size - 100ish
tstat     = 'F';            % 'F' for Granger's F-Test, 'chi2' for chi-squared
alpha     = 0.05;           % significance level
mhtc      = 'FDR';          % multiple hypothesis test correction (see routine 'significance')
k = 1:nobs;

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


tempX = X;
clear X
%% DO BIP STUFF HERE


% % lumped minus some other linear combo of unused channels
% coeffVec = [0.4; 0.3; 0.2; 0.1];
% multMat = repmat(coeffVec,[1 nobs nTrials]);
% ref = sum(multMat.*tempX(13:16,:,:));
% top = mean(tempX(1:4,:,:))-ref;
% mid = mean(tempX(5:8,:,:))-ref;
% dep = mean(tempX(9:12,:,:))-ref;

% % lumped minus some other linear combo, randomly selected
randoCoeffs = rand(4,1);
coeffVec = randoCoeffs/sum(randoCoeffs);
multMat = repmat(coeffVec,[1 nobs nTrials]);
ref = sum(multMat.*tempX(13:16,:,:));
top = mean(tempX(1:4,:,:))-ref;
mid = mean(tempX(5:8,:,:))-ref;
dep = mean(tempX(9:12,:,:))-ref;

% using channels actually in the analysis...
% top = mean(tempX(1:4,:,:))-tempX(11,:,:);
% mid = mean(tempX(5:8,:,:))-tempX(11,:,:);
% dep = mean(tempX(9:12,:,:))-tempX(11,:,:);

% bipolar subtraction
top = tempX(1,:,:)-tempX(3,:,:);
mid = tempX(5,:,:)-tempX(7,:,:);
dep = tempX(9,:,:)-tempX(11,:,:);


X = [top;mid;dep];
useBip = 0;

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

%% Model order estimation
momax = 30;
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

%% "Vertical" regression GC calculation
tic
modelOrder = 14;
shuffleCount = 10;
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

AVals = zeros([nvars, nvars, modelOrder, enobs]);
SIGVals = zeros([nvars, nvars, enobs]);

% loop through evaluation points
useVec = 35:55;
for e = useVec%1:enobs
    j = evalPoints(e);
    fprintf('window %d of %d at time = %d',e,enobs,j);

    [AT,SIG] = tsdata_to_var(X(:,j-wnobs+1:j,:),modelOrder,regmode);
    if isbad(AT)
        fprintf(2,' *** skipping - VAR estimation failed\n');
        continue
    end
    AVals(:,:,:,e) = AT;
    SIGVals(:,:,e) = SIG;

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
disp('Beginning Trial Shuffle runs')
clear i

dataPerms = zeros([shuffleCount,size(X)]); % scrambled input data
specGC_perm = zeros([shuffleCount,size(specGC)]); % scrambled GCs
timeGC_perm = zeros([shuffleCount, size(timeGC)]);
badCalcsPerm = zeros(shuffleCount,enobs);    % record where calculations fail


dispNullDists = ternaryOp(exist('dispNullDists','var'),dispNullDists,0); % default if needed
shuffleCount = ternaryOp(exist('shuffleCount','var'),shuffleCount,5); % default if needed

% trial shuffle and GC calcs

for k=1:shuffleCount
    % scramble the trials while retaining channel and temporal structure
    disp(['Shuffle ' num2str(k)])
    permX = zeros(size(X));
    for q = 1:size(X,1) % for each channel
        permX(q,:,:) = X(q,:,randperm(size(X,3)));  % dims: depth (channel), time, trial
    end
    dataPerms(k,:,:,:) = permX;

    
    % "Vertical" regression GC calculation
    specGCtemp = zeros(enobs,nvars,nvars,nBins);          % dims:  time, eq1, eq2, freq
    timeGCtemp = zeros(nvars, nvars, enobs);              % dims:  eq1, eq2, time
        
    % loop through evaluation points (where actual GC calculation is performed)
    for e = useVec
        j = evalPoints(e);
        fprintf('PERM %d: window %d of %d at time = %d',k,e,enobs,j);
        
        % convert time series data in window to VAR
        [A,SIG] = tsdata_to_var(permX(:,j-wnobs+1:j,:),modelOrder,regmode);
        if isbad(A)
            fprintf(2,' *** skipping - VAR estimation failed\n');
            badCalcsPerm(k,e) = 1;
            continue
        end
        
        % calculate autocovariance from VAR
        [G,info] = var_to_autocov(A,SIG);
        if info.error
            fprintf(2,' *** skipping - bad VAR (%s)\n',info.errmsg);
            badCalcsPerm(k,e) = 1;
            continue
        end
        if info.aclags < info.acminlags % warn if number of autocov lags is too small (not a show-stopper)
            fprintf(2,' *** WARNING: minimum %d lags required (decay factor = %e)',info.acminlags,realpow(info.rho,info.aclags));
            badCalcsPerm(k,e) = 1;
        end
        
        % Calculate GC for window
        a  = autocov_to_spwcgc(G,fres);
        specGCtemp(e,:,:,1:size(a,3)) = a;
        timeGCtemp(:,:,e) = squeeze(sum(specGCtemp(e,:,:,1:floor(nBins/5)),4))./(nBins/5); % up to ~60Hz
        if isbad(a,false)
            fprintf(2,' *** skipping - GC calculation failed\n');
            badCalcsPerm(k,e) = 1;
            continue
        end

        fprintf('\n');
    end
    
    specGC_perm(k,:,:,:,:) = specGCtemp;
    timeGC_perm(k,:,:,:) = timeGCtemp;

    toc
end

nullDistGC = squeeze(mean(specGC_perm))+2*squeeze(std(specGC_perm,0,1));
specGC = max(0,real(specGC));


%% plot pretty results

numVar = size(specGC,2);
maxGC = greatestMax(specGC(:,:,:,1:fres/2));

%SD Plot
figure()
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(specGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
            ylabel('Frequency (Hz)')
            axis xy
            axis([0 1 0 100])
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

%% frequency plot 

specGCShrunk = squeeze(sum(specGC,1))/length(useVec);
if length(size(nullDistGC)) > 3
    nullDistShrunk = squeeze(sum(nullDistGC,1))/length(useVec);
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

