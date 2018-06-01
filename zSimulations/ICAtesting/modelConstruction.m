%% Construct multiple non-stationary VAR time series
% model is derived from MVGC toolbox demo, var9model. Output is X
%
% Coefficients for the model are set in A, which is setup such that each
% "page" (along the 3rd dimension) of A can be multiplied with a vector of
% the channel's at the given page's lag to yield the AR equation. 
%
% For example, if we had a 2 variable 2nd order model defined by:
%
% A(:,:,1) = [0.5 0.2;
%            [0   -0.4];
% A(:,:,2) = [-0.3 0.1;
%             0    0.6];
%
% The model is produced by doing A(:,:,i)*[ X1[n-i] ; X2[n-i] ] for each
% page, i
%
% This would produce a model defined by the equations:
% X1[n] = 0.5*X1[n-1] + 0.2*X2[n-1] - 0.3*X1[n-2] + 0.1*X2[n-2]
% X2[n] = -0.4*X2[n-1] + 0.6*X2[n-2]
%
% WGN is then added to each equation based on the covaraince matrix, SIGT
% (which for simplicity of this demo, Barnett and Seth selected a 
% "minimal var" model using just an identity matrix for the covariance)
%
%

%% Parameters
tic
clear

ntrials   = 100;     % number of trials
dur       = 2;       % duration of trial (s)
fs        = 1000;    % sampling frequency (Hz)
seed      = 0;       % random seed (0 for unseeded) - 0

nobs  = dur*fs;      % number of observations per trial

%% VAR model construction
rng_seed(seed);

nvars = 9; % number of variables
p = 3;     % model order

A = zeros(nvars,nvars,p); % VAR coefficient matrix (dims: equation x affector x lag)

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


SIGT = eye(nvars); % covariance matrix; minimal VAR uses identity matrix

X = var_to_tsdata(A,SIGT,nobs,ntrials); % time series construction