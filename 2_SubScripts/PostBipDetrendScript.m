%% PostBipDetrendScript
% detrends data after bipolar subtraction

% extract and detrend data

% split up in and out targ data
intargTrials = ([data.inTarg] == 1);
outtargTrials = ([data.inTarg] == 0);
trialSetup = intargTrials; % so 1 = intarg, 0 = outtarg
intargData = data(intargTrials);
outtargData = data(outtargTrials);

% get size of matrices for detrend
tempTargDataIn = zeros([size(intargData(1).targlfpbip) length(intargData)]);
tempSaccDataIn = zeros([size(intargData(1).sacclfpbip) length(intargData)]);
tempTargDataOut = zeros([size(outtargData(1).targlfpbip) length(outtargData)]);
tempSaccDataOut = zeros([size(outtargData(1).sacclfpbip) length(outtargData)]);

% detrend based on in and out targ
disp('Detrending intarg')
tempTargDataIn = homeDetrend(reshape([intargData.targlfpbip],size(tempTargDataIn)));
tempSaccDataIn = homeDetrend(reshape([intargData.sacclfpbip],size(tempSaccDataIn)));
disp('Detrending outtarg')
tempTargDataOut = homeDetrend(reshape([outtargData.targlfpbip],size(tempTargDataOut)));
tempSaccDataOut = homeDetrend(reshape([outtargData.sacclfpbip],size(tempSaccDataOut)));

% reassign data
intargCount = 1;
outtargCount = 1;
for i=1:length(data)
    disp(['Reassigning trial ' num2str(i)])
    flag = trialSetup(i);
    if flag % intarg
    	data(i).targlfpmat = squeeze(tempTargDataIn(:,:,intargCount));
    	data(i).sacclfpmat = squeeze(tempSaccDataIn(:,:,intargCount));
    	intargCount = intargCount+1;
    else % outtarg
    	data(i).targlfpmat = squeeze(tempTargDataOut(:,:,outtargCount));
    	data(i).sacclfpmat = squeeze(tempSaccDataOut(:,:,outtargCount));
    	outtargCount = outtargCount+1;
    end
end


disp('Finished detrending')

