%% Trial shuffle runs.
% Run after performing GC normally!
%
% Establish GC null distribution by shuffling trials while retaining 
% temporal and spacial structure then finding GC
% For example, superficial -> mid GC normally is calculated by:
%
% trial 1 sup -> trial 1 mid
% trial 2 sup -> trial 2 mid
% ...
% trial 167 sup -> trial 167 mid
% 
% though this shuffle will now calculate GC by:
%
% trial 1 sup -> trial 24 mid
% trial 2 sup -> trial 152 mid
% etc.
% 

disp('Beginning Trial Shuffle runs')
clear i

dataPerms = zeros([shuffleCount,size(X)]); % scrambled input data
specGC_perm = zeros([shuffleCount,size(specGC)]); % scrambled GCs
timeGC_perm = zeros([shuffleCount, size(timeGC)]);
badCalcsPerm = zeros(shuffleCount,enobs);    % record where calculations fail


dispNullDists = ternaryOp(exist('dispNullDists','var'),dispNullDists,0); % default if needed
shuffleCount = ternaryOp(exist('shuffleCount','var'),shuffleCount,5); % default if needed

%% trial shuffle and GC calcs

for k=1:shuffleCount
    % scramble the trials while retaining channel and temporal structure
    disp(['Shuffle ' num2str(k)])
    permX = zeros(size(X));
    for q = 1:size(X,1) % for each channel
        permX(q,:,:) = X(q,:,randperm(size(X,3)));  % dims: depth (channel), time, trial
    end
    dataPerms(k,:,:,:) = permX;

    
    % "Vertical" regression GC calculation
    specGCtemp = nan(enobs,nvars,nvars,nBins);          % dims:  time, eq1, eq2, freq
    timeGCtemp = nan(nvars, nvars, enobs);              % dims:  eq1, eq2, time
        
    % loop through evaluation points (where actual GC calculation is performed)
    for e = 1:enobs
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
        
%         pval_t = mvgc_pval(timeGCtemp,modelOrder,wnobs,nTrials,1,1,nvars-2,tstat);
%         sig_t  = significance(pval_t,alpha,mhtc);
        
        fprintf('\n');
    end
    
    specGC_perm(k,:,:,:,:) = specGCtemp;
    timeGC_perm(k,:,:,:) = timeGCtemp;
    
    
    
    % Display results
    if dispNullDists
        % plot results
        numVar = size(specGC,2);
        maxGCtemp = greatestMax(specGCtemp);
        
        %SD Plot
        figure(68+3*k-2)
        for i=1:numVar
            for j=1:numVar
                if i~=j
                    subplot(numVar,numVar,(i-1)*numVar+j);
                    imagesc(specTime,freqs,squeeze(specGCtemp(:,i,j,:))', [0, maxGCtemp]) % why do I need to invert this? imagesc is weird :(
                    ylabel('Frequency (Hz)')
                    axis xy
                    axis([0 1 0 50])
                    %axis([-inf, inf, 0, 25])
                    colormap jet
                    set(gca, 'CLim', [0,maxGCtemp]);
                    hold on
                    plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
                    hold off
                end
            end
        end
        
        % for reference in figure
        subplot(numVar, numVar, numVar^2)
        set(gca, 'CLim', [0,maxGCtemp]);
        c = colorbar;
        c.Label.String = 'GC';
        colorbar
        title(['intarg = ' num2str(inTargVal)])
        xlabel(['Max GC = ' num2str(maxGCtemp) '      points/eval = ' num2str(pointsPerEval)])
        ylabel(['Model order = ' num2str(modelOrder)])
        
        %     %TD Plot
        %     maxTimeGC = greatestMax(timeGCtemp);
        %     figure(68+3*k-1)
        %     for i=1:numVar
        %         for j=1:numVar
        %             if i~=j
        %                 subplot(numVar,numVar,(i-1)*numVar+j);
        %                 plot(timeTime,squeeze(timeGCtemp(i,j,:)), 'LineWidth', 3)
        %                 ylabel('Granger Causality')
        %                 axis xy
        %                 axis([-inf, inf, 0, 1.2*maxTimeGC])
        %             end
        %         end
        %     end
        %
        %     % for reference in figure
        %     subplot(numVar, numVar, numVar^2)
        %     title(['intarg = ' num2str(inTargVal)])
        %     xlabel(['Max GC = ' num2str(maxGCtemp) '      points/eval = ' num2str(pointsPerEval)])
        %
    end
    toc
end
