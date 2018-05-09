%% Filtering script
% applies all desired filters to data
% takes rawDataGoodChan -> lfp

% Setup high pass filter frequency and coeffs
highpNorm = highpFreq/(rawFs/2);
[bHigh,aHigh] = butter(floor(totalHPFilterOrder/2),highpNorm,'high');

% Setup notch filter frequencies and coeffs
if harmonicNotch == 1
    notchfreqs = [notchFreq:notchFreq:5*notchFreq 321.7]/(rawFs/2);
elseif harmonicNotch == 2 % visually inspected noise bands
    notchfreqs = [60, 161, 300, 321.8, 483, 14276, 14919.6]/(rawFs/2);
else
    notchfreqs = [notchFreq 321.7]/(rawFs/2);
end
nnotch = length(notchfreqs);
anotch = zeros(nnotch,2*notchOrder+1); bnotch = zeros(nnotch,2*notchOrder+1);
for i=1:nnotch
   [bnotch(i,:),anotch(i,:)] = iirnotch(notchfreqs(i),notchfreqs(i)/30);
end

% Filter raw data for each trial
tic
for i=1:length(data)
    disp(['trial ' num2str(i)])
    rawData = double(data(i).rawDataGoodChan);

    notchedData = rawData;

    % Apply notch and high pass filters
    if doNotch
        for n=nnotch:-1:1
            notchedData = filtfilt(bnotch(n,:),anotch(n,:),notchedData')';
        end
    end
    
    if doHP
        data(i).lfp = filtfilt(bHigh,aHigh,notchedData')';
    else
        data(i).lfp = notchedData;
    end
    
end

disp(['notched @ ' num2str(notchfreqs*rawFs/2)])
