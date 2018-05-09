potato = zeros(3,1001,10);

t = 1:1001;

for i = 1:3
    for j = 1:10
        potato(i,:,j) = (j/i)*sin(2*pi*1/(rand*100)*t)+47;
    end
end

% demean by mvgc method
potato2 = demean(potato);

% demean by my method
potato3 = potato - repmat(mean(potato,3),1,1,size(potato,3));

% double demean (what's currently done)
potato4 = demean(potato3);


channelNum = randi(3);
trialNum = randi(10);

figure
subplot(4,1,1)
plot(squeeze(potato(channelNum, :, trialNum)))
hold on
plot(squeeze(mean(potato(channelNum,:,:),3)));
title('Original')
legend('Result','Mean Over Trials')

subplot(4,1,2)
plot(squeeze(potato2(channelNum, :, trialNum)))
hold on
plot(squeeze(mean(potato2(channelNum,:,:),3)));
title('MVGC demean')

subplot(4,1,3)
plot(squeeze(potato3(channelNum, :, trialNum)))
hold on
plot(squeeze(mean(potato3(channelNum,:,:),3)));
title('my Demean')

subplot(4,1,4)
plot(squeeze(potato4(channelNum, :, trialNum)))
hold on
plot(squeeze(mean(potato4(channelNum,:,:),3)));
title('both')

