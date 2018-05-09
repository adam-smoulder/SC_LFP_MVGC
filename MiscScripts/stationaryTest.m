%newX = diff(X,1,2);
newX = X;
%newX2 = X;
lags = 0;
startTimes = 0:10:800;


if length(lags) > 1 % if looking at diff lags, shows overall data
    results = nan([size(newX,1), size(newX,3), length(lags)]);
    pvals = nan(size(results));
    for i = 1:size(newX,1)
        for j = 1:size(newX,3)
            [results(i,j,:), pvals(i,j,:), ~,~,~] = kpsstest(newX(i,:,j),'lags',lags);
        end
    end
else % if only 1 lag, looks at individual time windows
    results = nan([size(newX,1), size(newX,3), length(startTimes)]);
    for i = 1:size(newX,1) % layer (sup/mid/deep)
        for j = 1:size(newX,3) % trial
            for k = 1:length(startTimes) % window in time
                results(i,j,k) = kpsstest(newX(i,startTimes(k)+1:startTimes(k)+200,j));
            end
        end
    end
end
resultsTotal = sum(sum(~results,2),1); % 0 means stationary from KPSS
disp(['Number of stationary: ' num2str(resultsTotal)])
disp(['Number of total obvs per lag/window: '...
    num2str(size(newX,3)*size(newX,1))])
