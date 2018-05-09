
% run fullBipMVGCPipeline for the fullBip data, both saccade onset and
% target onset
load('bb_sc_080415_fullBip.mat')
cueString = 'targlfp';  % cue to use
fullBipMVGCTemp

close all
clear


load('bb_sc_080415_fullBip.mat')
cueString = 'sacclfp';  % cue to use
fullBipMVGCTemp

close all
clear



% run timeMVGCPipeline on 083017 preproc data, just for saccade onset
load('ND_bl_sc_083017_preproc_det_2_HP_1_Notch_2_1000hz.mat')
timeMVGCPipeline

close all
clear




