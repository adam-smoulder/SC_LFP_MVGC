figs = [66];
numVars = 8;

isSgramsNotGC = 0;

clim = 0.225;
maxY = 100;


for q = figs
    figure(q)
    if ~isSgramsNotGC
        for i = 1:numVars
            for j = 1:numVars
                if numVars~=1
                    subplot(numVars,numVars,j+(i-1)*numVars)
                else
                    subplot(1,2,1)
                end
                axis xy
                axis([0 1 0 maxY])
                if clim
                    set(gca, 'CLim', [0,clim]);
                end
            end
        end
    else
        for i = 1:numVars
            if numVars~=1
                subplot(numVars,2,2*i-1)
            else
                subplot(1,2,1)
            end
            axis xy
            axis([0.2 1 0 maxY])
            if clim
                set(gca, 'CLim', [0,clim]);
            end
            
        end
    end
end
