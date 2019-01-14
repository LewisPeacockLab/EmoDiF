function [pcit_data] = emodif_prepare_pcit_data(target,sublist)

% Data format for P-CIT
% nTrials rows x 6 cols
% 1 - Subject ID
% 2 - Trial number
% 3 - Category
% 4 - Predictor variable
% 5 - Dependent variable
% 6 - Net effect cluster

%target = words, scenes, w-s, s-w 

%%bring