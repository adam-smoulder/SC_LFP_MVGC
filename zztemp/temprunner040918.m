% BAD channels by dataset (do NOT use these):
%
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

namesOfTheDay = [...
    'bb_sc_031315_mcell_spikelfp_cSC.mat';
    'bl_sc0723151_mcell_spikelfp_cSC.mat';
    'bb_sc_070915_mcell_spikelfp_cSC.mat';
    'bl_sc0723152_mcell_spikelfp_cSC.mat';
    'bb_sc_071215_mcell_spikelfp_cSC.mat';
    'bl_sc_031115_mcell_spikelfp_cSC.mat';
    'bb_sc_080415_mcell_spikelfp_cSC.mat';
    'bl_sc_071415_mcell_spikelfp_cSC.mat';
    'bb_sc_082917_mcell_spikelfp_cSC.mat';
    'bl_sc_083017_mcell_spikelfp_cSC.mat';
    'bb_sc_121415_mcell_spikelfp_cSC.mat';
    'bl_sc_112515_mcell_spikelfp_cSC.mat';
    ...
    ];

dasChannels = {...
    [1 3; 6 8; 10 12];
    [1 3; 5 7; 9 11];
    [1 3; 6 8; 10 12];
    [1 3; 5 7; 9 11];
    [1 3; 6 8; 10 12];
    [1 3; 6 8; 10 12];
    [1 3; 5 7; 9 11];
    [1 3; 6 8; 10 12];
    [1 3; 5 7; 9 11];
    [2 4; 6 8; 10 12];
    [1 3; 5 7; 9 11];
    [1 3; 5 7; 9 11];
    };


lNames = size(namesOfTheDay,1);

for theNameIWant = 1:lNames
    fname = namesOfTheDay(theNameIWant,:);
    channelsToUse = dasChannels{theNameIWant};
    disp(['loading ' fname]);
    preprocTemp
    clearvars -except namesOfTheDay lNames theNameIWant dasChannels fname
    disp(['FINISHED ' fname])
    close all
end

