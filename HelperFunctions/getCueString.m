function [ cueString ] = getCueString( cueType, signalType )
%getCueString returns the cueString from the given cueType and signalType
%   cueType: 'saccade', 'target', or 'go' - cue that data is centered
%   around in time
%   signalType: 'lfp' or 'spike' - may be open to combinations in future

switch cueType
    case 'saccade'
        cueString = strcat('sacc',signalType);
    case 'go'
        cueString = strcat('go',signalType);
    case 'target'
        cueString = strcat('targ',signalType);
    otherwise
        ME = MException('MyComponent:noSuchVariable','%s is not a cue name',cueType);
        throw(ME)
end
end

