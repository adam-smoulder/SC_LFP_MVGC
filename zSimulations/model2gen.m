% model 2 gen

t=linspace(0,dur,dur*fs);

%Generate the noise parameters
Var1=0.25;                      %Noise Variances
Var2=0.2;
Var3=0.30;
Var4=0.2;
%Initialize the output
N=length(t);
X = randn(N,4,ntrials);
%X = X-repmat(2*sin(2*pi*80*t)',[1 4 nTrials]); % add line noise if desired
X(:,1,:) = Var1*X(:,1,:);
X(:,2,:) = Var2*X(:,2,:);
X(:,3,:) = Var3*X(:,3,:);
X(:,4,:) = Var4*X(:,4,:);
X(2,1,:) = X(2,1,:) + 0.53*X(1,1,:);
X(2,2,:) = X(2,2,:) + 0.3*X(1,2,:);
X(2,3,:) = X(2,3,:) + 0.8*X(1,3,:);
X(2,4,:) = X(2,4,:) + 0.707*X(1,4,:);
for i = 3:N
   X(i,1,:) = X(i,1,:) + 0.5*X(i-1,1,:)-0.1*X(i-2,1,:)+0.16*X(i-2,3,:)+0.2*X(i-2,4,:);
   X(i,2,:) = X(i,2,:) + 0.5*X(i-1,2,:)-0.2*X(i-2,2,:)+0.1*X(i-1,3,:)+0.1*X(i-1,4,:);
   X(i,3,:) = X(i,3,:) + 0.8*X(i-1,3,:)-0.2*X(i-2,3,:);  %+0.5*X(i-2,2,:)
   X(i,4,:) = X(i,4,:) + .7*X(i-1,4,:)-0.5*X(i-2,4,:);
end

if makeSomeNoise == 1
    VarQ = mean(Var1, Var2, Var3, Var4);
    Q = noisePower*VarQ*repmat(randn(N,1),[1 4 ntrials]); % 0 mean side noise process
    X = X+Q;
elseif makeSomeNoise == 2
    R = zeros(size(X));
    maxptmag = greatestMax(X)+greatestMax(-X)-1;
    for i = 1:ntrials
        pt1 = maxptmag*(rand-1/2)+1;
        pt2 = maxptmag*(rand-1/2)+1;
        R(:,:,i) = noisePower*(squeeze(R(:,:,i))+...
            [repmat(linspace(pt1,pt2,N/2),[4 1])'; repmat(linspace(pt2,pt1,N/2),[4 1])']); % v shape
    end
    X = X+R;
end

% so Q is zero mean random white gaussian noise
% R is nonzero mean, linear, same for channels but changes every trial

%X = X-repmat(2*sin(2*pi*80*t)',[1 2 nTrials]); % add line noise if desired

tempX = zeros(nvars,nobs,ntrials);
X = reshape(X,size(tempX));

