% load data from trial shuffle, then run
% made for visualization of trial shuffle data beside normal lfps

figure(1)
trialnum = randi(size(X,3));
for i=1:3
    subplot(3,2,2*i-1)
    hold on
    plot(X(i,:,trialnum))
    %plot(X(i,:,1))
    plot(mean(X(i,:,:),3))
    axis([0,1000,-50,50])
end

figure(1)
for i = 1:3
    subplot(3,2,2*i)
    hold on
    plot(squeeze(dataPerms(1,i,:,trialnum)))
    plot(mean(squeeze(dataPerms(1,i,:,:)),2))
%     plot(squeeze(dataPerms(2,i,:,trialnum)))
%     plot(mean(squeeze(dataPerms(2,i,:,:)),2))
%     plot(squeeze(dataPerms(3,i,:,trialnum)))
%     plot(mean(squeeze(dataPerms(3,i,:,:)),2))
%     plot(squeeze(dataPerms(4,i,:,trialnum)))
%     plot(mean(squeeze(dataPerms(4,i,:,:)),2))
%     plot(squeeze(dataPerms(5,i,:,trialnum)))
%     plot(mean(squeeze(dataPerms(5,i,:,:)),2))
    axis([0,1000,-50,50])
end