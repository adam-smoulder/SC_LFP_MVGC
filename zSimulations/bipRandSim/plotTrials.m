myDataWhole = specGC_cont;

meanwave = squeeze(mean(myDataWhole));
medwave = squeeze(median(myDataWhole));

figure
for i = 1:3
    for j = 1:3
        if i ~=j
            myData = squeeze(myDataWhole(:,i,j,:));
            [coeff, score, eigenvals, tstat, percentExp, mu] = pca(squeeze(myData),'Centered',false);
            recon1 = mean(score(:,1)*coeff(:,1)',1);
            subplot(3,3,3*(i-1)+j)
            hold on
            for k = 1:20:1000
                plot(squeeze(myDataWhole(k,i,j,:)))
            end
            plot(squeeze(meanwave(i,j,:)),'r.-','Linewidth',4)
            %plot(recon1,'b.-','Linewidth',3)
            %plot(squeeze(medwave(i,j,:)),'g.-','Linewidth',2)
            axis([0,4000,0,1])
        end
    end
end
