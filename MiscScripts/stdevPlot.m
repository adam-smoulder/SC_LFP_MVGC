
%% plot results
numVar = size(stdSpecGC,2);
maxGC = greatestMax(stdSpecGC(:,:,:,1:fres/2));

%SD Plot
figure(66)
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(stdSpecGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
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

toc