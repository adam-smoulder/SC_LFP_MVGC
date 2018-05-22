

n = 9; % number of variables
p = 3; % model order

A = zeros(n,n,p);

A(:,:,1) = [
    0.0114         0    0.4088    0.4236         0         0         0         0         0;
   -0.3394   -0.0192         0    0.2799         0         0         0    0.3085         0;
         0         0    0.0194    0.3437         0         0         0         0         0;
         0    0.4302         0    0.0188         0         0         0         0         0;
         0         0    0.2851         0    0.0027         0    0.3160         0         0;
         0         0         0         0    0.3039   -0.0020         0         0         0;
         0         0         0         0         0    0.3030   -0.0186         0         0;
         0         0         0         0         0         0         0   -0.0084    0.3477;
         0         0         0         0         0         0    0.3037    0.2810   -0.0208
];

A(:,:,2) = [
    0.0148         0    0.2590    0.1965         0         0         0         0         0;
    0.2232   -0.0070         0    0.1779         0         0         0   -0.2383         0;
         0         0   -0.0058    0.2008         0         0         0         0         0;
         0    0.2103         0    0.0012         0         0         0         0         0;
         0         0    0.1597         0    0.0065         0    0.1989         0         0;
         0         0         0         0    0.2062    0.0177         0         0         0;
         0         0         0         0         0    0.1895   -0.0008         0         0;
         0         0         0         0         0         0         0   -0.0032    0.1381;
         0         0         0         0         0         0    0.1947   -0.1718    0.0068
];

A(:,:,3) = [
    0.0076         0    0.1434    0.1787         0         0         0         0         0;
    0.1082   -0.0065         0    0.1351         0         0         0    0.2021         0;
         0         0    0.0050    0.1826         0         0         0         0         0;
         0    0.2090         0   -0.0037         0         0         0         0         0;
         0         0    0.1943         0    0.0097         0   -0.0698         0         0;
         0         0         0         0    0.1347   -0.0081         0         0         0;
         0         0         0         0         0    0.2270   -0.0004         0         0;
         0         0         0         0         0         0         0   -0.0014    0.1397;
         0         0         0         0         0         0    0.2020    0.1417   -0.0022
];

AT = A;

nvars = size(AT,1); % number of variables

SIGT = eye(nvars);

X = var_to_tsdata(AT,SIGT,nobs,ntrials);

xMax = greatestMax(X);
xMin = -greatestMax(-X);

if makeSomeNoise == 1
    base = randn(N,ntrials);
    VarQ = (xMax-xMin)/greatestMax(base);
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
    R = permute(R,[3 1 2]);
    X = X+R;
end
