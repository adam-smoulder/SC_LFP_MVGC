%% Realign script
% aligns lfp data about target onset and saccade onset
% takes lfp -> sacclfpmat and targlfpmat

for i=1:length(data)
    disp(['trial ' num2str(i)])
    trialData = double(data(i).lfp);
        
    % get data for cue times w.r.t. start of trial
    targCode = data(i).targCode;
    goCode = data(i).goCode;
    measCode = data(i).measCode;
    srt = data(i).srts;
    currstates = data(i).stateTransitions;
    
    targtime = double(currstates(2,currstates(1,:)==targCode))*30;
    gotime = double(currstates(2,currstates(1,:)==goCode))*30;
    sacctime = (double(currstates(2,currstates(1,:)==measCode))+srt)*30;
    
    % realign data for target onset
    targtimevec = (data(i).targtimevec(1)*30):(data(i).targtimevec(2)*30);
    data(i).targlfpmat = trialData(:,floor(targtime+targtimevec));
    
    % realign data for go cue
    if(~isempty(gotime))
        gotimevec = data(i).gotimevec(1)*30:data(i).gotimevec(2)*30;
        data(i).golfpmat = trialData(:,floor(gotime+gotimevec));
    end
    
    % realign data for saccade onset
    sacctimevec = data(i).sacctimevec(1)*30:data(i).sacctimevec(2)*30;
    data(i).sacclfpmat = trialData(:,floor(sacctime+sacctimevec));
end


