% Load 16ChanGC data, then select desired channels for affected
% (vars2checki) and affector (vars2checkj).
%% plot results
vars2checki=9:16;
vars2checkj=9:16;
numVar = length(vars2checki);
maxGC = greatestMax(specGC(:,:,:,1:fres/2));

%SD Plot
figure(66)
for i=vars2checki
    for j=vars2checkj
        if i~=j
            subplot(numVar,numVar,(mod(i-1,8)+1-1)*numVar+mod(j-1,8)+1);
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


% %TD Plot
% maxTimeGC = greatestMax(timeGC);
% figure(67)
% for i=1:numVar
%     for j=1:numVar
%         if i~=j            
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             plot(timeTime,squeeze(timeGC(i,j,:)), 'LineWidth', 3)
%             ylabel('Granger Causality')
%             axis xy
%             axis([-inf, inf, 0, 1.2*maxTimeGC])
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar, numVar, numVar^2)
% title(['intarg = ' num2str(inTargVal)])
% xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])

toc
