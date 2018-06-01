rows = 1:12;
cols = 1:12;
numVar = length(rows);
%%
%SD Plot
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
title(['row: ' num2str(min(rows)) '-' num2str(max(rows))])
ylabel(['col: ' num2str(min(cols)) '-' num2str(max(cols))])

% for reference in figure
subplot(numVar,numVar, (numVar^2-floor(length(rows))/2))
xlabel('Frequency')
subplot(numVar,numVar,(numVar*floor((numVar-1)/2)+1))
ylabel('GC Power')
subplot(numVar,numVar,2)
legend('Actual','Null')
subplot(numVar, numVar, numVar^2)
xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
ylabel(['Model order = ' num2str(modelOrder)])

%% establish the significance image and plot it
crunchedGC = squeeze(sum(specGCShrunk,3))/500;
crunchedNull = squeeze(sum(nullDistShrunk,3))/500;
sigImage = (crunchedGC > crunchedNull);

subplot(1,3,1)
imagesc(crunchedGC)
set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
colorbar
subplot(1,3,2)
imagesc(crunchedNull)
set(gca, 'CLim', [0,greatestMax(crunchedGC)]);   
colorbar
subplot(1,3,3)
imagesc(sigImage)


