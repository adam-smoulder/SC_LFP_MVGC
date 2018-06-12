clearvars -except specGC_groundTruth specGC_noPollCAR specGC_randBipMean specGC_CAR specGC_bip shit newOrder trueGC

fileToUse = 'beforeTestClose.mat';
load(fileToUse)
% 
% if makeSomeNoise == 1
%     base = randn(nobs,ntrials);
%     VarQ = noisePower*(xMax-xMin)/(greatestMax(base)+greatestMax(-base));
%     baseRepeat = permute(repmat(base,[1 1 nvars]), [3 1 2]);
%     Q = VarQ*baseRepeat-mean([xMax xMin])+1; % 0 mean (?) white noise process
%     trueX = trueX+Q;
% end
trueX = trueX_NN;
modelOrder = 3; % or 13 for ding model maybe?
nobs = size(trueX,2);
ntrials = size(trueX,3);
fres = 4000;
fs = 200;
regmode = 'OLS';
windowSize = nobs-modelOrder;
pointsPerEval = windowSize;

dnobs     = 0;       % initial observations to discard per trial - 0
regmode   = 'OLS';   % VAR model estimation regression mode ('OLS', 'LWR' or empty for default - 'OLS' ?
nlags     = 50;       % number of lags in VAR model

dur = nobs/fs;
tnobs = nobs+dnobs;                 % total observations per trial for time series generation
k = 1:tnobs;                        % vector of time steps
t = (k-dnobs-1)/fs;                 % vector of times
%% RANDOM BIP STUFF SETUP HERE
%makeSomeNoise = 1;
%noisePower = 10;
%theSetup = eye(9);

CAR = 0; % 0 = no CAR, 1 = normal CAR, 2 = no-pollution CAR
bipRand = 0; % 0 = no bip, 1 = biprand, 2 = bip normal

noPoll = 1; % no pollution in biprand 4
randStyle = 1; % 1 = diffSet, 2 = diffTrial, 3 = diffTime, 4 = diffChan
indMod = 0; % mod for diffChan
randPick = ''; % 'normal' for normally distributed about mean (randbip)
fullOrNot = 'full'; % for rand bip
shuffleCount = 5;
numruns = 100;
% modelOrder = 3;

useSets = 0; % to ignore subtraction of specific sets for CAR
sets = [1 2; 1 2; 3 4; 3 4; 5 6; 5 6; 7 8; 7 8; 8 9;];

if CAR==1 % normal
    tempX = [trueX(:,:,:)];
    trueX = trueX-repmat(mean(tempX),[9,1,1]);
elseif CAR==2 % no-pollution
    coeffOrder = 1:9;%randperm(size(trueX,1));%
    carVal = zeros(size(trueX)); % the effective reference to subtract
    for i = 1:size(trueX,1)
        if useSets
            coeffs = 1/(size(trueX,1)-2)*ones(size(trueX,1),1);
            coeffs(squeeze(sets(i,:))) = 0;
        else
            coeffs = 1/(size(trueX,1)-1)*ones(size(trueX,1),1);
            coeffs(coeffOrder(i)) = 0;
        end
        carVal(i,:,:) = squeeze(sum(trueX.*coeffs));
    end
    trueX = trueX-carVal;
end
% CAR does best; both upweight and downweight make it worse
% Increasing upweight removes existing GC
% increasing downweight introduces more spurious GC
% Nonpolluted shows...exactly the same? is that expected?

if bipRand == 1 % bip rand
    CVs = [];
    for datCount = 1:numruns
        clear X
%         X = zeros(size(trueX)-[13,0,0]);
%         refChans = 13:16;
%         topChans = 1:4;
%         midChans = 5:8;
%         deepChans = 9:12;
        X = zeros(size(trueX)-[6,0,0]);
        if strcmp(fullOrNot,'full')
            refChans = 1:9;
        elseif strcmp(fullOrNot,'not')
            refChans = 7:9;
        elseif strcmp(fullOrNot,'oddball')
            refChans = [1:9 1:6 1:6 1:6];
        end
        topChans = 1:2;
        midChans = 3:4;
        deepChans = 5:6;
        
        nref = length(refChans);
        
        % make reference for removal using random weights (same over trials)
        if randStyle == 1 % diff datasets
            nref = length(refChans);
            if strcmp(randPick,'normal')
                randoCoeffs = (1/nref)+(1/nref)/4*randn(length(refChans),1);
                %randoCoeffs = randoCoeffs.*[1;1;1;1;1;1;2;2;2];
            else
                randoCoeffs = rand(length(refChans),1);
            end
            % randoCoeffs = (theSetup(datCount,:))'; % fake; use for trying stuff
            coeffVec = randoCoeffs/sum(randoCoeffs);
            CVs = [CVs; coeffVec'];
            multMat = repmat(coeffVec,[1 nobs]);
        end
        
        % for each trial, randbip it
        for i = 1:ntrials
            disp(['Bip replacement ' num2str(i) ' of ' num2str(ntrials)])
            if randStyle == 1 % diff datasets
                ref = squeeze(sum(multMat.*squeeze(trueX(refChans,:,i)))); %ONLY FOR FROZEN CVs over trial
            elseif randStyle == 2 % diff trials
%             % make reference for removal using random weights, diff trials
                    nref = length(refChans);
                    if strcmp(randPick,'normal')
                        randoCoeffs = (1/nref)+(1/nref)/4*randn(length(refChans),1);
                    else
                        randoCoeffs = rand(length(refChans),1);
                    end                %randoCoeffs = (theSetup(datCount,:))'; % fake; use for trying stuff
                coeffVec = randoCoeffs/sum(randoCoeffs);
                %CVs = [CVs; coeffVec'];
                multMat = repmat(coeffVec,[1 nobs]);
                ref = squeeze(sum(multMat.*trueX(refChans,:,i)));
            elseif randStyle == 3 % diff time
                CVforTime = [];
                for j = 1:nobs
                    nref = length(refChans);
                    if strcmp(randPick,'normal')
                        randoCoeffs = (1/nref)+(1/nref)/4*randn(length(refChans),1);
                    else
                        randoCoeffs = rand(length(refChans),1);
                    end
                    coeffVec = randoCoeffs/sum(randoCoeffs);
                    %CVs = [CVs; coeffVec'];
                    CVforTime = [CVforTime coeffVec];
                end
                multMat = CVforTime;
                ref = squeeze(sum(multMat.*trueX(refChans,:,i)));
            elseif randStyle == 4 % time and channel
                ref = zeros(nref,nobs);
                for j = 1:nobs
                    if strcmp(randPick,'normal')
                        randoCoeffs = (1/nref)+(1/nref)/4*randn(nref, nref);
                    else
                        randoCoeffs = rand(nref,nref);
                    end
                    if noPoll
                        randoCoeffs = (1-eye(nref)).*randoCoeffs;
                    end
                    coeffMat = randoCoeffs./sum(randoCoeffs);
                    %CVs = [CVs; coeffVec'];
                    %CVforTime = [CVforTime coeffMat];
                    
                    for k = 1:nref
                        if indMod
                            randChanToSub = randi(nref);
                            while randChanToSub == k
                                randChanToSub = randi(nref);
                            end
                            ref(k,j) = trueX(refChans(randChanToSub),j,i);
                        else
                            ref(k,j) = squeeze(sum(coeffMat(:,k).*trueX(refChans,j,i)));
                        end
                    end
                end
                
                % subtract the reference and reassign
                topTarg = squeeze(mean(trueX(topChans,:,i)-ref(topChans,:)));
                midTarg = squeeze(mean(trueX(midChans,:,i)-ref(midChans,:)));
                deepTarg = squeeze(mean(trueX(deepChans,:,i)-ref(deepChans,:)));
                X(:,:,i) = [topTarg; midTarg; deepTarg];
            end
            
            if randStyle ~= 4 % method 4 is hackish atm...
                % subtract the reference and reassign
                topTarg = squeeze(mean(trueX(topChans,:,i)))-ref;
                midTarg = squeeze(mean(trueX(midChans,:,i)))-ref;
                deepTarg = squeeze(mean(trueX(deepChans,:,i)))-ref;
                X(:,:,i) = [topTarg; midTarg; deepTarg];
            end
        end
                
        % DONT FORGET TO DETREND
        X = homeDetrend(X);
        
        
        
        
        % "Vertical" regression GC calculation
        nvars = size(X,1);
        
        wnobs = modelOrder+windowSize;                  % number of observations in "vertical slice"
        evalPoints    = wnobs : pointsPerEval : nobs;   % GC evaluation points
        enobs = length(evalPoints);                     % number of GC evaluations
        f1 = zeros(enobs,nvars,nvars,fres);             % (time, eq1, eq2, freq)
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
    specGCstd = squeeze(std(specGC_cont,1,1));
    specGCvar = specGCstd.^2;
    nullDistGC = squeeze(mean(nullDistGC_cont));
    
else
    if bipRand == 2 % bip sub
%         top = trueX(1,:,:)-trueX(3,:,:);
%         mid = trueX(5,:,:)-trueX(7,:,:);
%         deep = trueX(9,:,:)-trueX(11,:,:);
          top = trueX(1,:,:)-trueX(2,:,:);
          mid = trueX(3,:,:)-trueX(4,:,:);
          deep = trueX(5,:,:)-trueX(6,:,:);
    else % just lumping if biprand isn't 1 or 2
        top = mean(trueX(1:2,:,:));
        mid = mean(trueX(3:4,:,:));
        deep = mean(trueX(5:6,:,:));
    end
    
    X = [top;mid;deep];
    
    
    X = homeDetrend(X);

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
    nullDistGC = squeeze(mean(specGC_perm))+4*squeeze(std(specGC_perm,0,1));
    
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
maxGC = greatestMax(trueGC);

rows = 1:numVar;
cols = 1:numVar;

errors = sum(abs(trueGC-squeeze(specGC)),3);
spurError = sum(errors([3 7 8]));
lostError = sum(errors([2 4 6]));

figure()
for i=1:length(rows)
    for j=1:length(cols)
        if rows(i)~=cols(j)
            subplot(numVar,numVar,(i-1)*numVar+j);
            plot(freqs,squeeze(specGCShrunk(rows(i),cols(j),:)),'Linewidth',2)
            hold on
            plot(freqs,squeeze(nullDistShrunk(rows(i),cols(j),:)),'LineWidth',2)
            plot(freqs,squeeze(trueGC(rows(i),cols(j),:)),'Linewidth',1)
            title(['error = ' num2str(errors(i,j))])
            axis([0 fs/2 0 maxGC])
        end
    end
end
subplot(numVar,numVar,1)
%title(['row: ' num2str(min(rows)) '-' num2str(max(rows)) ' col: ' num2str(min(cols)) '-' num2str(max(cols))])
xlabel('Frequency')

% for reference in figure
subplot(numVar,numVar,2)
legend('Predicted','Null','Truth')
subplot(numVar, numVar, numVar^2)
title(['SpurErr = ' num2str(spurError) ', lostErr = ' num2str(lostError)])
xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
ylabel(['Model order = ' num2str(modelOrder)])


%%
% establish the significance image and plot it
% crunchedGC = squeeze(sum(specGCShrunk,3))/(fs/2);
% crunchedNull = squeeze(sum(nullDistShrunk,3))/(fs/2);
% sigImage = (crunchedGC > crunchedNull);
% 
% figure
% subplot(1,3,1)
% imagesc(crunchedGC)
% title('GC')
% set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
% colorbar
% subplot(1,3,2)
% imagesc(crunchedNull)
% title('Null Dist')
% set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
% colorbar
% subplot(1,3,3)
% imagesc(sigImage)
% title('GCs above null dist')


%%
% <mvgc_demo_nonstationary.html back to top>