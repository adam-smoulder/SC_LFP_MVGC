%% script for testing effectiveness of ICA for removal of added reference
%  see modelConstruction script for how the model was made; you can rerun
%  it yourself if you have the MVGC toolbox added to your path!
%  This script assumes fastICA is in MATLAB's path

tic
clear
%% Parameters

load('modelForICATest.mat') % modelConsruction
% data for analysis is in X, dims: channel x time x trial (realization)

addNoise     = 2;    % 0 = no noise, 1 = WGN, 2 = linear non-zero mean 2-step
noisePower   = 1;    % magnitude of noise added w.r.t. input signal (SNR)
icaVarThresh = inf;  % maximum variance between weights for the selected IC

%% Adding noise if desired to X 
xMax = squeeze(max(max(max((X)))));
xMin = squeeze(min(min(min((X)))));

if addNoise == 1    % zero mean WGN
    unscaledRef = randn(nobs,ntrials);                   % zero-mean WGN, std of 1
    signalStd = std(reshape(X,[1 numel(X)]));            % stdev of signal's distribution
    scaledRef = signalStd*unscaledRef;                   % so SNR is 1
    R = permute(repmat(scaledRef,[1 1 nvars]), [3 1 2]); % reshape and copy over channels
elseif addNoise == 2 % non-zero mean 2-step linear reference
    R = zeros([nobs nvars ntrials]);
    maxptmag = xMax-xMin-0.5;
    for i = 1:ntrials
        pt1 = maxptmag*(rand-1/2)+1; % random start/stop point
        pt2 = maxptmag*(rand-1/2)+1; % random midpoint
        R(:,:,i) = squeeze(R(:,:,i))+...
            [repmat(linspace(pt1,pt2,nobs/2),[nvars 1])'; ...
            repmat(linspace(pt2,pt1,nobs/2),[nvars 1])']; % v shape b/w pt1 -> pt2 -> pt1
    end
    R = permute(R,[2 1 3]);
else
    R = zeros(size(X)); 
end
% plot(R(1,:,1)) % see an example of the noise

X = X+noisePower*R;


%% perform ICA
varDims = cell([1, ntrials]); % variance between channel weights for ICs (just for reference)
predRef = zeros(size(X));  % predicted reference signal with IC and A
%predRef = zeros([nobs ntrials]); % predicted with weight and original
% I think these ^^ two methods should be similar? I think another step is
% needed for the commented method, as it's off by a lot in its DC component

for i = 1:ntrials % run ICA separately on each trial
    Xi = squeeze(X(:,:,i));
    [S, A, W] = fastica(Xi); % perform ICA for the trial
    [varA, minVarDim] = min(var(A));          % find dimension with smallest variance
    if varA < icaVarThresh                    % if under thresh, remove suspected ref.
        predRef(:,:,i) = A(:,minVarDim)*S(minVarDim,:);
        %predRef(:,i) =  W(minVarDim,:)*Xi;
    end
    disp(['Trial ' num2str(i) ' of ' num2str(ntrials) ' complete'])
    varDims(i) = {var(A)}; % saves variance spread of A for later reference
end

%X = X-removedRef; % if desired, subtract predicted reference from data


timeToComplete = toc;

%% plotting reconstruction vs. actual reference
% for visualization, I'll plot the reconstruction of the predicted
% reference for 3 channels along with the actual reference signal

trialsToUse = [randi(ntrials) randi(ntrials) randi(ntrials)];
chansToUse = [1 4 8];
legendString = 'legend(';

t = linspace(0,dur,nobs*length(trialsToUse));
refForPlot = reshape(squeeze(R(1,:,trialsToUse)),size(t));
predRefForPlot = reshape(squeeze(predRef(:,:,trialsToUse)),[nvars length(t)]);
%predRefForPlot = reshape(squeeze(predRef(:,trialsToUse)),size(t));

figure
hold on
for i = chansToUse
    plot(t,predRefForPlot(i,:))
    legendString = [legendString '"chan' num2str(i) '",'];
end
plot(t,refForPlot, 'k--', 'LineWidth',2)
title(['Trials ' num2str(trialsToUse)])
xlabel('Time (s)')
ylabel('Voltage (V)')
legendString = [legendString(1:end) ' "R")'];
eval(legendString)
