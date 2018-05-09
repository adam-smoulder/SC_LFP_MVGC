% ADAM BL SACC NEEDS TO BE RUN WITH NO DETREND MO62!

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

namesOfTheDay2 = [...
% 'ND_bb_sc_031315_preproc_det_0_HP_1_Notch_2_1000hz.mat';
% 'ND_bb_sc_070915_preproc_det_0_HP_1_Notch_2_1000hz.mat';
% 'ND_bb_sc_071215_preproc_det_0_HP_1_Notch_2_1000hz.mat';
% 'ND_bb_sc_080415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
% 'ND_bb_sc_082917_preproc_det_0_HP_1_Notch_2_1000hz.mat';
% 'ND_bb_sc_121415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc0723151_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc0723152_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc_031115_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc_071415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc_083017_preproc_det_0_HP_1_Notch_2_1000hz.mat';
'ND_bl_sc_112515_preproc_det_0_HP_1_Notch_2_1000hz.mat';
    ...
    ];


lNames = size(namesOfTheDay2,1);

% for theNameIWant2 = 1:lNames
%     fname = squeeze(namesOfTheDay2(theNameIWant2,:));
%     load(fname)
%     cueString = 'targlfp';
%     disp(['loading ' fname]);
%     mvgcTemp2
%     clearvars -except namesOfTheDay2 lNames theNameIWant2 fname
%     disp(['FINISHED ' fname])
%     close all
% end


for theNameIWant2 = 1:lNames
    fname = squeeze(namesOfTheDay2(theNameIWant2,:));
    load(fname)
    cueString = 'sacclfp';
    disp(['loading ' fname]);
    mvgcTemp2
    clearvars -except namesOfTheDay2 lNames theNameIWant2 fname
    disp(['FINISHED ' fname])
    close all
end




