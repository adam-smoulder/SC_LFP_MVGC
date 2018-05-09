function [ outputCSD ] = lfp2csd( inputData )
% converts lfps to csd using the iCSD2D toolbox
%   

ntrials = size(inputData, 3);
nt = size(inputData,2);
nx = 2;
ny = size(inputData,1);
dx = 0.1;
dy = 0.1;
VX = 1:dx:nx;
VY = 1:dy:ny;
fileName = ['csdRaw' num2str(randi(1000))];
method = 'lin';

allPotentials = zeros([ntrials nt nx ny]); % dims: trial, time, x, y
for i=1:ny
    allPotentials(:,:,1,i) = squeeze(inputData(i,:,:))';
    allPotentials(:,:,2,i) = squeeze(inputData(i,:,:))';
end

initlin2d(fileName, nx, ny, dx, dy, 1);
csd = zeros([ntrials, nt, nx, ny]);
out = zeros([ntrials, nt, length(VX), length(VY)]);
load(['rawcsd_data/' fileName])
for i = 1:ntrials
    csd(i,:,:,:) = icsd2d(squeeze(allPotentials(i,:,:,:)), F);
    out(i,:,:,:) = interp2d(squeeze(csd(i,:,:,:)), VX, VY, method);
end

outputCSD = squeeze(out(:,:,1,:));

end

