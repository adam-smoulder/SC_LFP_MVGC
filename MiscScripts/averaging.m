% script for averaging together spectral GC sets

normie = 0; % if 1, normalizes max GC of each dataset to 1 pre-averaging

% names = [...
%     'bb_sc_070915_targlfp_inTarg_1_MO_65_wind_100_1000hz.mat';
%     'bb_sc_071215_targlfp_inTarg_1_MO_65_wind_100_1000hz.mat';
%     'bb_sc_080415_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723151_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723152_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc_071415_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat'];

% names = [...
%         'bb_sc_031315_targlfp_MO_62_intarg_1.mat';
%     'bb_sc_071215_targlfp_inTarg_1_MO_65_wind_100_1000hz.mat';
%     'bb_sc_080415_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723151_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723152_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc_071415_targlfp_inTarg_1_MO_59_wind_100_1000hz.mat'];

% names = [...
%     'bb_sc_070915_sacclfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     %'bb_sc_071215_sacclfp_inTarg_1_MO_65_wind_100_1000hz.mat';
%     'bb_sc_080415_sacclfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723151_sacclfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc0723152_sacclfp_inTarg_1_MO_59_wind_100_1000hz.mat';
%     'bl_sc_071415_sacclfp_inTarg_1_MO_59_wind_100_1000hz.mat';];

thenamesForHere = [...
%     'bb_sc_031315_sacclfp_MO_62_intarg_1.mat'; % saccade
%     'bb_sc_070915_sacclfp_MO_62_intarg_1.mat';
%     'bb_sc_071215_sacclfp_MO_62_intarg_1.mat';
%     'bb_sc_080415_sacclfp_MO_62_intarg_1.mat';
%     'bb_sc_082917_sacclfp_MO_62_intarg_1.mat';
%     'bb_sc_121415_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc_071415_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc0723151_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc0723152_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc_031115_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc_112515_sacclfp_MO_62_intarg_1.mat';
%     'bl_sc_083017_sacclfp_MO_62_intarg_1.mat';

    'ND_bb_sc_031315_targlfp_MO_25_intarg_1.mat'; % some acausal target?
    'ND_bb_sc_070915_targlfp_MO_25_intarg_1.mat';
    'ND_bb_sc_071215_targlfp_MO_25_intarg_1.mat';
    'ND_bb_sc_080415_targlfp_MO_25_intarg_1.mat';
    'ND_bb_sc_082917_targlfp_MO_25_intarg_1.mat';
    'ND_bb_sc_121415_targlfp_MO_25_intarg_1.mat';
    'ND_bl_sc0723151_targlfp_MO_25_intarg_1.mat';
    'ND_bl_sc0723152_targlfp_MO_25_intarg_1.mat'; % virtually no sacc response?
    'ND_bl_sc_031115_targlfp_MO_25_intarg_1.mat';
    'ND_bl_sc_071415_targlfp_MO_25_intarg_1.mat';
    'ND_bl_sc_083017_targlfp_MO_25_intarg_1.mat'; % noise band around 180?
%     %'ND_bl_sc_112515_targlfp_MO_62_intarg_1.mat'; % acausal target???
% % %     
%     'ND_bb_sc_031315_sacclfp_MO_25_intarg_1.mat';
%     'ND_bb_sc_070915_sacclfp_MO_25_intarg_1.mat';
%     'ND_bb_sc_071215_sacclfp_MO_25_intarg_1.mat';
%     'ND_bb_sc_080415_sacclfp_MO_25_intarg_1.mat';
%     'ND_bb_sc_082917_sacclfp_MO_25_intarg_1.mat';
    %'ND_bb_sc_121415_sacclfp_MO_25_intarg_1.mat';
%     'ND_bl_sc0723151_sacclfp_MO_25_intarg_1.mat';
%     %'ND_bl_sc0723152_sacclfp_MO_62_intarg_1.mat';
%     'ND_bl_sc_031115_sacclfp_MO_25_intarg_1.mat';
%     'ND_bl_sc_071415_sacclfp_MO_25_intarg_1.mat';
%     'ND_bl_sc_083017_sacclfp_MO_25_intarg_1.mat';
%     'ND_bl_sc_112515_sacclfp_MO_25_intarg_1.mat';
    ...
%     'bb_sc_031315_targlfp_MO_62_intarg_1.mat'; % target
%     'bb_sc_070915_targlfp_MO_62_intarg_1.mat';
%     'bb_sc_071215_targlfp_MO_62_intarg_1.mat';
%     'bb_sc_080415_targlfp_MO_62_intarg_1.mat';
%     'bb_sc_082917_targlfp_MO_62_intarg_1.mat';
%     'bb_sc_121415_targlfp_MO_62_intarg_1.mat';
%     'bl_sc_071415_targlfp_MO_62_intarg_1.mat';
%     'bl_sc0723151_targlfp_MO_62_intarg_1.mat';
%     'bl_sc0723152_targlfp_MO_62_intarg_1.mat';
%     'bl_sc_031115_targlfp_MO_62_intarg_1.mat';
%     'bl_sc_112515_targlfp_MO_62_intarg_1.mat';
%     'bl_sc_083017_targlfp_MO_62_intarg_1.mat'...
%       'DET0_bb_sc_031315_targlfp_MO_62_1hz.mat'; % detrend zero
%       'DET0_bb_sc_070915_sacclfp_MO_62_1hz.mat';
%       'DET0_bb_sc_071215_sacclfp_MO_62_1hz.mat';
%       'DET0_bb_sc_080415_sacclfp_MO_62_1hz.mat';
%       'DET0_bb_sc_082917_sacclfp_MO_62_1hz.mat';
%       'DET0_bb_sc_121415_targlfp_MO_62_1hz.mat';
    ];
 
specGCVals = cell(size(thenamesForHere,1),1);
gcSizes = zeros(1,size(thenamesForHere,1));

for counter=1:size(thenamesForHere,1)
    load(thenamesForHere(counter,:));
    gcSizes(counter) = size(specGC,1);
    specGC(~(specGC >= 0)) = 0; % convert negatives and nans to 0
    specGCVals{counter} = specGC;
    if counter~=size(thenamesForHere,1)
        clear specGC
    end
    disp(['Loaded ' num2str(counter)])
end

sizeToUse = max(gcSizes);
gcMatrix = zeros([size(thenamesForHere,1) sizeToUse, size(specGC,2), size(specGC,3), size(specGC,4)]);
clear specGC


for counter=1:size(thenamesForHere,1)
    specGC = specGCVals{counter};
    if size(specGC,1) ~= sizeToUse % if wrong size, concatenate blank
        difference = sizeToUse-size(specGC,1);
        specGCnew = zeros(size(specGC)+[difference,0,0,0]);
        specGCnew(1:end-difference,:,:,:) = specGC;
        specGCnew(end-difference+1:end,:,:,:) = specGC(end-difference+1,:,:,:);
        specGC = specGCnew;
    end
    if normie
        specGC = specGC/greatestMax(specGC);
    end
    gcMatrix(counter,:,:,:,:) = specGC;
    
    clear specGC
end


avgSpecGC = squeeze(mean(gcMatrix));
expAvgSpecGC = squeeze(mean(exp(gcMatrix),1));
sterrorSpecGC = squeeze(std(gcMatrix,1))/sqrt(size(specGCVals,1));
expStdSpecGC = squeeze(std(exp(gcMatrix),1));
medianSpecGC = squeeze(median(gcMatrix));

specGC = avgSpecGC;

%% plot results
numVar = size(avgSpecGC,2);
maxGC = greatestMax(avgSpecGC(:,:,:,1:fres/2));

%mean Plot
figure(1)
for counter=1:numVar
    for j=1:numVar
        if counter~=j
            subplot(numVar,numVar,(counter-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(avgSpecGC(:,counter,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
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

%% plot std error
numVar = size(sterrorSpecGC,2);
maxGC = greatestMax(sterrorSpecGC(:,:,:,1:fres/2));

%SD Plot
figure(2)
for counter=1:numVar
    for j=1:numVar
        if counter~=j
            subplot(numVar,numVar,(counter-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(sterrorSpecGC(:,counter,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
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


% %% plot exp mean
% numVar = size(expAvgSpecGC,2);
% maxGC = greatestMax(expAvgSpecGC(:,:,:,1:fres/2));
% minGC = min(min(min(min(expAvgSpecGC(:,:,:,1:fres/2)))));
% 
% %exp mean Plot
% figure(3)
% for i=1:numVar
%     for j=1:numVar
%         if i~=j
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             imagesc(specTime,freqs,squeeze(expAvgSpecGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
%             ylabel('Frequency (Hz)')
%             axis xy
%             axis([0 1 0 50])
%             %axis([-inf, inf, 0, 25])
%             colormap jet
%             set(gca, 'CLim', [minGC,maxGC]);
%             hold on
%             plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
%             hold off            
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar, numVar, numVar^2)
% set(gca, 'CLim', [0,maxGC]);
% c = colorbar;
% c.Label.String = 'GC';
% colorbar
% title(['intarg = ' num2str(inTargVal)])
% xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
% ylabel(['Model order = ' num2str(modelOrder)])
% 
% toc
% 
% %% plot exp stdev
% numVar = size(expStdSpecGC,2);
% maxGC = greatestMax(expStdSpecGC(:,:,:,1:fres/2));
% minGC = min(min(min(min(expStdSpecGC(:,:,:,1:fres/2)))));
% 
% 
% %geometric SD Plot
% figure(4)
% for i=1:numVar
%     for j=1:numVar
%         if i~=j
%             subplot(numVar,numVar,(i-1)*numVar+j);
%             imagesc(specTime,freqs,squeeze(expStdSpecGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
%             ylabel('Frequency (Hz)')
%             axis xy
%             axis([0 1 0 50])
%             %axis([-inf, inf, 0, 25])
%             colormap jet
%             set(gca, 'CLim', [minGC,maxGC]);
%             hold on
%             plot(startTime*dur*ones(1,100),1:100,'r--','LineWidth',2)
%             hold off            
%         end
%     end
% end
% 
% % for reference in figure
% subplot(numVar, numVar, numVar^2)
% set(gca, 'CLim', [0,maxGC]);
% c = colorbar;
% c.Label.String = 'GC';
% colorbar
% title(['intarg = ' num2str(inTargVal)])
% xlabel(['Max GC = ' num2str(maxGC) '      points/eval = ' num2str(pointsPerEval)])
% ylabel(['Model order = ' num2str(modelOrder)])
% 
% toc

%% plot median
numVar = size(medianSpecGC,2);
maxGC = greatestMax(medianSpecGC(:,:,:,1:fres/2));

%median Plot
figure(5)
for counter=1:numVar
    for j=1:numVar
        if counter~=j
            subplot(numVar,numVar,(counter-1)*numVar+j);
            imagesc(specTime,freqs,squeeze(medianSpecGC(:,counter,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
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
