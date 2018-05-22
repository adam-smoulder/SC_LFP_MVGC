% model 3 generation 
% should show 1 -> 2,3,4 and 4<->5

N = length(t);
n = 5; % number of variables
p = 3; % max model order

r = sqrt(2);

A = zeros(n,n,p);

%coeff syntax: (a,b,c)
% a = target eq, b = origination of signal eq, c = index difference (-1,-2, -3)

A(1,1,1) =  0.95*r; %1->self
A(1,1,2) = -0.9025; %1->self

A(2,1,2) =  0.5; %1->2
A(3,1,3) = -0.4; %1->3

A(4,1,2) = -0.5; %1->4
A(4,4,1) =  0.25*r; %4->self
A(4,5,1) =  0.25*r; %5->4

A(5,4,1) = -0.25*r; %4->5
A(5,5,1) =  0.25*r; %5->self

AT = A;

nvars = size(AT,1); % number of variables

SIGT = eye(nvars);

X = var_to_tsdata(AT,SIGT,nobs,ntrials);
xMax = greatestMax(X);
xMin = -greatestMax(-X);

if makeSomeNoise == 1
    base = randn(N,ntrials);
    VarQ = noisePower*(xMax-xMin)/(greatestMax(base)+greatestMax(-base));
    baseRepeat = permute(repmat(base,[1 1 nvars]), [3 1 2]);
    Q = VarQ*baseRepeat-mean([xMax xMin])+1; % 0 mean side noise process
    X = X+Q;
elseif makeSomeNoise == 2
    R = zeros([nobs nvars ntrials]);
    maxptmag = noisePower*(greatestMax(X)+greatestMax(-X)-1);
    for i = 1:ntrials
        pt1 = maxptmag*(rand-1/2)+1;
        pt2 = maxptmag*(rand-1/2)+1;
        R(:,:,i) = squeeze(R(:,:,i))+...
            [repmat(linspace(pt1,pt2,N/2),[nvars 1])'; repmat(linspace(pt2,pt1,N/2),[nvars 1])']; % v shape
    end
    R = permute(R,[3 1 2]);
    X = X+R;
end
