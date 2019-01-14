function [p] = PCIT_parms()
%parameters for PCIT runs - this is for PCIT data preparation - this is
%changed EVERY TIME you want to generate new PCIT data


%% Analysis ID
% Give this a unique name that will be used to save the bootstrap data


%% TR window configuration - do we even need 
% these?  
% p.trStarts = 1:7;
% p.trLens   = 3:7;
% p.nBins    = 3;

%%% other parameters DO NOT CHANGE
p.nSubjects = 24;
p.nCats = 2;
p.nCol_pcit = 6;
p.nTrials = 60;
p.nTrials_SS = p.nSubjects * p.nTrials;


%category, and trial selection
% mask = 1 (VTC), 2 (HIP)

p.act  = '1';

switch p.act
    case '1'
        p.target = 'word';
    case '2'
        p.target = 'scene';
    case '3'
        p.target = 'face';
    case '4'
        p.target = 'object';
    case '5'
        p.target = 'rest';
    case '6'
        p.target = 'word-scene';
end
        


%% Bootstrap configuration
%
p.bsSampSubj = false;
p.nBoots     = 1;
p.nSubjects = 24;

%script directories
%on Astoria
p.scriptdir = '/Users/tw24955/EmoDIF/PCIT';

end