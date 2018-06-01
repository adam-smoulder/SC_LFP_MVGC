n = 16; % number of variables
p = 14; % model order

N = nobs;

load('chan16A41andSIG41.mat')
AT = A41;
nvars = size(AT,1); % number of variables
SIGT = SIG41;
% SIGT = eye(nvars); % minimal var

X = var_to_tsdata(AT,SIGT,nobs,ntrials);

xMax = greatestMax(X);
xMin = -greatestMax(-X);

if makeSomeNoise == 1
    base = randn(N,ntrials);
    VarQ = noisePower*(xMax-xMin)/greatestMax(base);
    baseRepeat = permute(repmat(base,[1 1 nvars]), [3 1 2]);
    Q = VarQ*baseRepeat-mean([xMax xMin])+1; % 0 mean side noise process
    X = X+Q;
elseif makeSomeNoise == 2
    R = zeros([nobs nvars ntrials]);
    maxptmag = greatestMax(X)+greatestMax(-X)-0.5;
    for i = 1:ntrials
        pt1 = maxptmag*(rand-1/2)+1;
        pt2 = maxptmag*(rand-1/2)+1;
        R(:,:,i) = squeeze(R(:,:,i))+...
            [repmat(linspace(pt1,pt2,N/2),[nvars 1])'; repmat(linspace(pt2,pt1,N/2),[nvars 1])']; % v shape
    end
    R = permute(R,[2 1 3]);
    X = X+noisePower*R;
end
