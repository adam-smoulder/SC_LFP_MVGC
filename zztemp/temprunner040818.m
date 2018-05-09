namesOfTheDay = [...
    'bb_sc_031315_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bb_sc_070915_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bb_sc_071215_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bb_sc_080415_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bb_sc_082917_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bb_sc_121415_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc_071415_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc0723151_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc0723152_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc_031115_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc_112515_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    'bl_sc_083017_preproc_det_2_HP_1_Notch_2_1000hz.mat';
    ...
%    'BipPerm_bb_sc_080415_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    'BipPerm_bb_sc_082917_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    'BipPerm_bb_sc_121415_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    'BipPerm_bl_sc0723151_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    'BipPerm_bl_sc0723152_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    'BipPerm_bl_sc_112515_filted_det_2_HP_1_Notch_2_1000hz.mat';
%    ...
%     'bb_sc_031315_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bb_sc_070915_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bb_sc_071215_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bb_sc_080415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bb_sc_082917_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bb_sc_121415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc_071415_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc0723151_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc0723152_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc_031115_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc_112515_preproc_det_0_HP_1_Notch_2_1000hz.mat';
%     'bl_sc_083017_preproc_det_0_HP_1_Notch_2_1000hz.mat'...
    ];


lNames = size(namesOfTheDay,1);

for theNameIWant = 1:lNames
    fname = squeeze(namesOfTheDay(theNameIWant,:));
    load(fname);
    cueString = 'targlfp';
    mvgcTemp
    clearvars -except namesOfTheDay lNames theNameIWant
    close all
    
%     load(names(thename,:));
%     cueString = 'sacclfp';
%     mvgcTemp
%     clearvars -except names lNames thename
    disp(['FINISHED PART ' num2str(theNameIWant)])
%     close all
end


