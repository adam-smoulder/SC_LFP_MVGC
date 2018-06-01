exper = reshape(removedRef, nvars, ntrials*nobs); % placing tiles horizontally
actual = reshape(R, nvars, ntrials*nobs); % placing tiles horizontally

for i = 1:nvars
    subplot(nvars,1,i)
    plot(exper(i,:))
    hold on
    plot(actual(i,:),'LineWidth',1)
    axis([30000 40000 -inf inf])
end