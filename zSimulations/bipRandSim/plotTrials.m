myDataWhole = squeeze(specGC_cont);
%load('NN_noBip')
specGC = trueGC;

meanwave = squeeze(mean(myDataWhole));
medwave = squeeze(median(myDataWhole));

figure
for i = 1:3
    for j = 1:3
        if i ~=j
            myData = squeeze(myDataWhole(:,i,j,:));
            [coeff, score, eigenvals, tstat, percentExp, mu] = pca(squeeze(myData));
            recon1 = mean(score(:,1)*coeff(:,1)',1);
            subplot(3,3,3*(i-1)+j)
            hold on
            for k = 1:2:100
                plot(squeeze(myDataWhole(k,i,j,:)))
            end
            plot(squeeze(meanwave(i,j,:)),'r.-','Linewidth',3)
            plot(squeeze(specGC(i,j,:)),'b.-','Linewidth',3)
            %plot(recon1,'b.-','Linewidth',3)
            %plot(squeeze(medwave(i,j,:)),'g.-','Linewidth',2)
            axis([0,4000,0,0.6])
        end
    end
end

%% plot spectra
spectra_NN = abs(fft(trueX_NN,[],2));
spectra_Q1 = abs(fft(trueX_Q1,[],2));
spectra_Q10 = abs(fft(trueX_Q10,[],2));
spectra_justQ1 = abs(fft(trueX_Q1-trueX_NN,[],2));

useThisOne = spectra_Q1;


meanwave = squeeze(mean(useThisOne,3));
medwave = squeeze(median(useThisOne,3));

nChan = size(meanwave,1);
freqs = linspace(0,fs/2,ceil(length(meanwave)/2));
figure
for i = 1:nChan
    subplot(nChan,1,i)
    plot(freqs,meanwave(i,1:end/2))
    hold on
    plot(freqs,medwave(i,1:end/2))
    if i == 1
        legend('mean','median')
    elseif i == 5
        ylabel('Amplitude')
    elseif i == 9
        xlabel('Frequency (Hz)')
    end
    axis([1 fs/2 0 max(meanwave(i,2:end/2))])
end




