% ADAM: load data appropriately, then individually run sections

%% Intarg confidence interval plotting
%load('my intarg data here')
inTarg_meanTimeGC_perm = squeeze(mean(timeGC_perm));
inTarg_stdTimeGC_perm = squeeze(std(timeGC_perm));
inTarg_meanSpecGC_perm = squeeze(mean(specGC_perm));
inTarg_stdSpecGC_perm = squeeze(std(specGC_perm));

numVar = nvars;

maxTimeGC = greatestMax(timeGC);

figure(66)
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            hold on
            plot(timeTime,squeeze(timeGC(i,j,:)), 'LineWidth', 3)
            plot(timeTime, squeeze(inTarg_meanTimeGC_perm(i,j,:)),'k', 'LineWidth',2)
            plot(timeTime, squeeze(inTarg_meanTimeGC_perm(i,j,:)+2.*inTarg_stdTimeGC_perm(i,j,:)), 'r--','LineWidth',2)
            plot(timeTime, squeeze(inTarg_meanTimeGC_perm(i,j,:)-2.*inTarg_stdTimeGC_perm(i,j,:)), 'r--','LineWidth',2)
            ylabel('Granger Causality')
            axis xy
            axis([-inf, inf, 0, 1.2*maxTimeGC])
        end
    end
end

% for reference in figure
subplot(numVar, numVar, numVar^2)
%title(['intarg = ' num2str(inTargVal)])
title('intarg = 1')
xlabel(['points/eval = ' num2str(pointsPerEval)])


%% Outtarg confidence interval
%load('my outtarg data here')
outTarg_meanTimeGC_perm = squeeze(mean(timeGC_perm));
outTarg_stdTimeGC_perm = squeeze(std(timeGC_perm));
outTarg_meanSpecGC_perm = squeeze(mean(specGC_perm));
outTarg_stdSpecGC_perm = squeeze(std(specGC_perm));

maxTimeGC = greatestMax(timeGC);
numVar = nvars;

figure(67)
for i=1:numVar
    for j=1:numVar
        if i~=j            
            subplot(numVar,numVar,(i-1)*numVar+j);
            hold on
            plot(timeTime,squeeze(timeGC(i,j,:)), 'LineWidth', 3)
            plot(timeTime, squeeze(outTarg_meanTimeGC_perm(i,j,:)),'k', 'LineWidth',2)
            plot(timeTime, squeeze(outTarg_meanTimeGC_perm(i,j,:)+2.*outTarg_stdTimeGC_perm(i,j,:)), 'r--','LineWidth',2)
            plot(timeTime, squeeze(outTarg_meanTimeGC_perm(i,j,:)-2.*outTarg_stdTimeGC_perm(i,j,:)), 'r--','LineWidth',2)
            ylabel('Granger Causality')
            axis xy
            axis([-inf, inf, 0, 1.2*maxTimeGC])
        end
    end
end

% for reference in figure
subplot(numVar, numVar, numVar^2)
%title(['intarg = ' num2str(inTargVal)])
title('intarg = 0')
xlabel(['points/eval = ' num2str(pointsPerEval)])