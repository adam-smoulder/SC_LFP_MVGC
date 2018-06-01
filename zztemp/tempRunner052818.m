% bb 070915: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bb 071215: 11, 14 -> use [1 3; 6 8; 10 12]
% bb 080415: none - > [1 3; 5 7; 9 11] 
% bb 031315:  11 14 (use [1 3; 5 7; 9 12])
% bb 121415: none [1 3; 5 7; 9 11]
% bb 082917: none [1 3; 5 7; 9 11]
% 24bb 012317: none
% 24bb 081617: none

% bl 071415: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 0723151: none ->  using [1 3; 5 7; 9 11]
% bl 0723152: none -> using [1 3; 5 7; 9 11]
% bl 031115: 2, 7, 11, 14 -> use [1 3; 6 8; 10 12]
% bl 112515:  none [1 3; 5 7; 9 11]
% bl 083017:  1 -> use [2 4; 6 8; 10 12]
% 24bl 071717: none
% 24bl 071718: none
% 24bl 071719: none

numruns = 10;

targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);

for datCount = 1:numruns
    fname = 'bb_sc_070915_mcell_spikelfp_cSC.mat';
    topChans = [1 3 4];
    midChans = [5 6 8];
    deepChans = [9 10 12];
    refChans = [13 15 16];
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns




targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bb_sc_071215_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = [9 10 12];
    refChans = [13 15 16];
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns

%% done up to here

targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bb_sc_080415_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns




targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bb_sc_031315_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = [9 10 12];
    refChans = [13 15 16];
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns





targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bb_sc_082917_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns



 

targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc_071415_mcell_spikelfp_cSC.mat';
    topChans = [1 3 4];
    midChans = [5 6 8];
    deepChans = [9 10 12];
    refChans = [13 15 16];
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns



targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc0723151_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns




targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc0723152_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns



targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc_031115_mcell_spikelfp_cSC.mat';
    topChans = [1 3 4];
    midChans = [5 6 8];
    deepChans = [9 10 12];
    refChans = [13 15 16];
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns




targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc_083017_mcell_spikelfp_cSC.mat';
    topChans = 2:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns



% these ones tend to take longer, sometimes fail...

targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bb_sc_121415_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns



targSpecGC_holder = zeros([89, 3, 3, 8001]);
saccSpecGC_holder = zeros([89, 3, 3, 8001]);
targNullGC_holder = zeros([89, 3, 3, 8001]);
saccNullGC_holder = zeros([89, 3, 3, 8001]);
for datCount = 1:numruns
    fname = 'bl_sc_112515_mcell_spikelfp_cSC.mat';
    topChans = 1:4;
    midChans = 5:8;
    deepChans = 9:12;
    refChans = 13:16;
    tempPreproc052818
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'targlfp';
    tempMVGC052818
    targSpecGC_holder = targSpecGC_holder+specGC;
    targNullGC_holder = targNullGC_holder+nullDistSpecGC;
    close all
    
    load(strcat(outputFileName,'.mat'));
    cueString = 'sacclfp';
    saccSpecGC_holder = saccSpecGC_holder+specGC;
    saccNullGC_holder = saccNullGC_holder+nullDistSpecGC;
    tempMVGC052818
    close all
end
clear specGC nullDistSpecGC
targSpecGC = targSpecGC_holder/numruns;
targNullGC = targNullGC_holder/numruns;
saccSpecGC = saccSpecGC_holder/numruns;
saccNullGC = saccNullGC_holder/numruns;
outputFileName3 = strcat('randBipMEAN_',fname(1:end-23),...
    '_MO_',num2str(modelOrder),...
    '_intarg_',num2str(inTargVal),'.mat');
save(outputFileName3,'-v7.3')

clearvars -except numruns

