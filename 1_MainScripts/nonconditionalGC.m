% % manually load data (currently 031315)
% % remove goodLFP = 0
% % sort out inTarg = 1 and inTarg = 0
% %
%
%
% % vars
% fs = 1000;                                  % sampling frequency
% window = fs/5;                               % size of window for FFT
% slide = fs/100;
% nFFT = 1024*4;                                % number of points for FFT
%
% % data input for target timepoint
% % for trial = 1:length(data)
% %     trialData = data(i);
% %     for i=1:length(trialData.channelOrder)
% %         sortedLFPdata(trial, trialData.channelOrder(i),:) = trialData.targlfpmat(i,:);
% %     end
% % end
%
%
% %%

%% GC Prep variables that may change between runs

cueString = 'sacclfp';  % cue to use
inTargVal = 1;          % intarg vs outtarg vs both

% Setup for GC

icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

% Parameters that should not change between runs
dnobs     = 0;          % initial observations to discard per trial - 0; we can ignore as needed
dur       = 1;          % duration of trial (s) - 1s

nobs  = dur*fs+1;       % number of observations in a trial
tnobs = nobs+dnobs;     % total observations per trial for time series generation
k = 1:tnobs;            % vector of time steps

monkey_data = fname(1:12);          %
axisLimit = 0;                      % if 0, axisLimit = maxGC (only used for subtractorExtractor)
downsampleFactor = 1000/fs;         % takes 1khz -> 250hz % HACK fix later


%% Isolate desired data

% inTarg = 0 -> out of target, inTarg = 1 -> in target, otherwise -> all
switch inTargVal
    case 0
        dataToUse = data([data.inTarg] == 0 & [data.goodLFP] == 1);
    case 1
        dataToUse = data([data.inTarg] == 1  & [data.goodLFP] == 1);
    otherwise
        dataToUse = data;
end

% specifically isolate the cue we want to analyze
[nChannels, nTime] = size(dataToUse(1).sacclfpmat);
nTrials = length(dataToUse);  % updated now for in vs. outtarg
X = nan(nChannels, nTime, nTrials);
if strcmp(cueString,'targlfp')
    X = reshape([dataToUse.targlfpmat],size(X));
    startTime = 0.4;
elseif strcmp(cueString,'sacclfp')
    X = reshape([dataToUse.sacclfpmat],size(X));
    startTime = 0.6;
else
    disp('Something went wrong...Check your cueString')
end

% remove the rest of the stuff from memory, we shouldn't need it!
clear data;
clear dataToUse;
%%

% isolate data desired to be analyzed
%dataToUse = sortedLFPdata(:,1:2:13,:);
%dataForGC = reshape(X,size(sortedLFPdata));

% load data, run through first two blocks of mvgcPipeline
oldX = X;
window = fs/5;                               % size of window for FFT
slide = fs/100;
nFFT = 1024*4;

dataForGC = permute(oldX,[2 1 3]); % rearrange for this analysis

%Generate the Spectral Density Matrix
tic;
Results = grangerCoordinator(dataForGC,fs,0,window,slide,nFFT);
toc

tempSpecGC = Results.GrangerCausality;
time = Results.time;
freqs = Results.frequency;

specGC = permute(tempSpecGC,[4 1 2 3]); % arrange back

%% plotting
numVar = size(specGC,2);
maxGC = greatestMax(specGC(:,:,:,1:length(freqs)/2));
pointsPerEval = slide;

%SD Plot
figure(66)
for i=1:numVar
    for j=1:numVar
        if i~=j
            subplot(numVar,numVar,(i-1)*numVar+j);
            imagesc(time,freqs,squeeze(specGC(:,i,j,:))', [0, maxGC]) % why do I need to invert this? imagesc is weird :(
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

% %% Display the Results
% figure(1); %Pairwise GC
% figure(2); %Coherence
% numVar = size(dataForGC,2);
% k = 1;
% for r = 1:numVar
%     for c = 1:numVar
%         if r == c
% %             figure(1);
% %             subplot(5,5,k);
% %             plot(Results.frequency,squeeze(Results.SpectralDensity(r,c,:)))
% %             set(gca,'ylim',[0 1])
% %             figure(2);
% %             subplot(5,5,k);
% %             plot(Results.frequency,squeeze(Results.SpectralDensity(r,c,:)))
% %             set(gca,'ylim',[0 1])
%         else
%             figure(1);
%             subplot(numVar,numVar,k);
%             plot(Results.frequency,squeeze(Results.GrangerCausality(r,c,:)))
%             set(gca,'ylim',[0 1])
%             figure(2);
%             subplot(numVar,numVar,k);
%             plot(Results.frequency,squeeze(Results.Coherence(r,c,:)))
%             set(gca,'ylim',[0 1])
%         end
%         k = k+1;
%     end
% end

%% save
disp('Saving...')

clear X oldX dataForGC tempSpecGC

outputFileName = strcat('GC_',fname(1:end-23),...
    '_',cueString,...
    '_inTarg_',num2str(inTargVal),...
    '_nonpar',...
    '_',num2str(fs),'hz.mat');

save(outputFileName,'-v7.3')

disp('MVGC Complete!')
