tempsup = squeeze(simplifiedCsdData(1,:,1));
tempmid = squeeze(simplifiedCsdData(1,:,2));
tempdeep = squeeze(simplifiedCsdData(1,:,3));

figure(1)
plot(tempsup)
hold on
plot(tempmid)
plot(tempdeep)
legend('sup','mid','deep')

trialnumpick = 43;

trial1csd = zeros(1001,151);
for i = 1:151
    trial1csd(:,i) = squeeze(csdData(trialnumpick,:,i));
end

trial1lfp = squeeze(sortedData(:,:,trialnumpick))';



% 
% figure(2)
% subplot(3,1,1)
% sum1 = zeros(1,1001);
% hold on
% for i = 11:21
%     plot(trial1(:,i));
%     sum1 = sum1+squeeze(trial1(:,i));
% end
% plot(sum1./11,'r--','LineWidth',2)
% hold off
% 
% figure(2)
% subplot(3,1,2)
% sum2 = zeros(1,1001);
% hold on
% for i = 51:61
%     plot(trial1(:,i));    
%     sum2 = sum2+squeeze(trial1(:,i));
% end
% plot(sum2./11,'r--','LineWidth',2)
% hold off
% 
% figure(2)
% subplot(3,1,3)
% sum3 = zeros(1,1001);
% hold on
% for i = 91:101
%     plot(trial1(:,i));
%     sum3 = sum3+squeeze(trial1(:,i));
% end
% plot(sum3./11,'r--','LineWidth',2)
% hold off
% 
% 
% %axis adjustment
% figure(2)
% for i=1:3
%     subplot(3,1,i)
%     axis([-inf, inf, -80000, 100000])
% end


figure(3)
for i = 1:4
    subplot(3,2,1)
    hold on
    plot(trial1csd(:,(i-1)*10+1));
    title('CSD')
    ylabel('chan 1-4')
    hold off
    
    subplot(3,2,2)
    hold on
    plot(trial1lfp(:,i))
    title('LFP')
    ylabel('chan 1-4')
    hold off
end

figure(3)
for i = 5:8
    subplot(3,2,3)
    hold on
    plot(trial1csd(:,(i-1)*10+1));  
    ylabel('chan 5-8')
    hold off
   
    subplot(3,2,4)
    hold on
    plot(trial1lfp(:,i))
    ylabel('chan 5-8')
    hold off
end

figure(3)
for i = 9:12
    subplot(3,2,5)
    hold on
    plot(trial1csd(:,(i-1)*10+1));
    ylabel('chan 9-12')
    hold off
    
    subplot(3,2,6)
    hold on
    plot(trial1lfp(:,i))
    ylabel('chan 9-12')
    hold off
end

%axis adjustment
figure(3)
for i=1:3
    subplot(3,2,(2*i)-1)
    axis([-inf, inf, -80000, 100000])
    subplot(3,2,(2*i))
    axis([-inf, 1001, -100, 100])
end


figure(4)
subplot(1,2,1)
hold on
imagesc(trial1csd(:,1:121)')
colormap jet
axis xy
axis([0,1001,1,121])
title('CSD')
ylabel('depth')
hold off

subplot(1,2,2)
hold on
imagesc(trial1lfp(:,1:13)')
colormap jet
axis xy
axis([0,1001,1,13])
title('LFP')
ylabel('depth (channel)')
hold off
