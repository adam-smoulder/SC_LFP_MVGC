% for quick use

potato = reshape(trialLfp',[1 numel(trialLfp)]);
w = linspace(-rawFs/2,rawFs/2,length(potato));
specPotato = abs(fftshift(fft(potato)));
plot(w,specPotato)
axis([0 inf -inf inf])

%%
w(find(specPotato > 40000))