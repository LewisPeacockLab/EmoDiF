function emodif_onset_update(subjNum)

%this is to update the SM of preview to disregard cues and to add 2TRs to
%instructions so that we can look at the effect of instructions.  this is
%to look at the effect of instruction and also to model with the onset and
%instruction (with onset being non-critical)

subjID = sprintf('emodif_%s',num2str(subjNum));
script_dir = pwd;

%read in old DFencode hi conident onsets
load(sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RF_SM_hiconf_onsets_%s.mat',subjID,subjID));

onset_instr{1,1} = onsets{1,1}+2;
onset_instr{1,2} = onsets{1,2}+2;
onset_instr{1,3} = onsets{1,3}+2;
onset_instr{1,4} = onsets{1,4}+2;

onset_item = onsets;
onsets = onset_instr;


outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RF_SM_hiconf_instr_onsets_%s.mat', subjID, subjID);
save(outfname,'names','onsets','durations');

names = {'for_SF_item','for_SF_instr','for_SR_item','for_SR_instr','rem_SF_item', 'rem_SF_instr', 'rem_SR_item', 'rem_SR_instr'};

onsets{1,1} = onset_item{1,1};
onsets{1,2} = onset_instr{1,1};
onsets{1,3} = onset_item{1,2};
onsets{1,4} = onset_instr{1,2};
onsets{1,5} = onset_item{1,3};
onsets{1,6} = onset_instr{1,3};
onsets{1,7} = onset_item{1,4};
onsets{1,8} = onset_instr{1,4};

durations = ones(length(onsets),1)';
durations= num2cell(durations);

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RF_SM_hiconf_combined_onsets_%s.mat', subjID, subjID);
save(outfname,'names','onsets','durations');



%emotional
load(sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RF_SM_hiconf_onsets_%s.mat',subjID,subjID));
onset_instr{1,1} = onsets{1,1}+2;
onset_instr{1,2} = onsets{1,2}+2;
onset_instr{1,3} = onsets{1,3}+2;
onset_instr{1,4} = onsets{1,4}+2;
onset_item = onsets;
onsets = onset_instr;

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RFemo_instr_onsets_%s.mat', subjID, subjID);
save(outfname,'names','onsets','durations');


names = {'for_neg_item','for_neg_instr','for_neu_item','for_neu_instr','rem_neg_item', 'rem_neg_instr','rem_neu_item','rem_neu_instr'};

onsets{1,1} = onset_item{1,1};
onsets{1,2} = onset_instr{1,1};
onsets{1,3} = onset_item{1,2};
onsets{1,4} = onset_instr{1,2};
onsets{1,5} = onset_item{1,3};
onsets{1,6} = onset_instr{1,3};
onsets{1,7} = onset_item{1,4};
onsets{1,8} = onset_instr{1,4};

durations = ones(length(onsets),1)';
durations= num2cell(durations);

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RFemo_combined_onsets_%s.mat', subjID, subjID);
save(outfname,'names','onsets','durations');
delete(sprintf('~/emodif_data/%s/behav/EmoDif_SPM_DFEncode_RFemo_hiconf_combined_onsets_%s.mat',subjID, subjID));

%read in old onsets
load(sprintf('~/emodif_data/%s/behav/EmoDif_mvpa_allregs.mat', subjID));

%error in 105-124 where for emotion 1 is -1
mvpa_regs.preview.emo(mvpa_regs.preview.emo<0) = 1;
mvpa_regs.DFEncode.emo(mvpa_regs.DFEncode.emo<0) = 1;
%%%NOTE THIS DOES NOT FIX THIS IN THE REGRESSORS - this needs to be put in
%%%all MVPA parse files too %%%%


%%%%% these variables ALL change with experiment %%%%
%write names and durations
if subjID == 101 | subjID == 102 | subjID == 103
    localizer_run_TR = 426;
else
localizer_run_TR = 342; %426 in subject 1, 2 and 3.  %342 for everyone else.  
preview_run_TR = 366;
DFencode_run_TR = 426;
end

names_local = mvpa_regs.localizer.cat_names; % easy %face scene object word rest, for first three subjects - word is emotional and neutral
names_DFE = mvpa_regs.DFEncode.cat_names;
names_preview = mvpa_regs.preview.cat_names;

% duration is 0 if we are modeling a delta - this goes up in number for
% boxcar models etc.  1 duration for each condition.
% durations = zeros(length(names),1)'; % easy **** NOW FOUND UNDER PHASES



%%%%%%%% BUILDING our onsets %%%%%%%%%%%%%%
% for each of these conditions:
%face scene object word rest
%%% Start Massive Indexing Effort or SMIE %%%
loc_trials = mvpa_regs.localizer.trial;
DFE_trials = mvpa_regs.DFEncode.trial;
preview_trials = mvpa_regs.preview.trial;


DFE_resp = mvpa_regs.DFEncode.subresp;
preview_resp = mvpa_regs.preview.subresp;

% 
% face_index = mvpa_regs.localizer.cat==1;
% scene_index = mvpa_regs.localizer.cat==2;
% object_index = mvpa_regs.localizer.cat ==3;

% if subjID == 101 || subjID == 102 || subjID == 103
% neg_word_index = mvpa_regs.localizer.cat==4;
% pos_word_index = mvpa_regs.localizer.cat ==5;   
% rest_index = mvpa.regs.localizer.cat ==6;
% else
%     
% word_index = mvpa_regs.localizer.cat==4;
% rest_index = mvpa_regs.localizer.cat ==5; 
% end

%onsets are all expanded - we need to find the beginning of each miniblock
%for the localizer 

%for exceptions - identify by trial number


%general parameters
betweenrun_rest = 3;
%localizer parameters
loc_miniblock_TR = 9;
loc_postblock_rest = 5;
%DFencode parameters

DFEncode_trialTR=7; %how many TRS are in one trial for DFencoding
preview_trialTR=6; %how many TRs are in one trial for preview
DFEncode_trialN = 60;
preview_trialN = 60;

%preview parameters

 trial_start_run1 = 1:preview_trialTR:(preview_trialTR*preview_trialN)/2;
    trial_start_run2 = (preview_trialTR*(preview_trialN/2)+betweenrun_rest)+1:preview_trialTR:(preview_trialTR*preview_trialN)+betweenrun_rest;
    
    %here are exceptions for preview
    if subjNum == 107 %
        bigger = find(trial_start_run1 > 162);
        for x = 1:length(bigger)
            trial_start_run1(bigger(x)) = trial_start_run1(bigger(x))+1;
        end
    elseif subjNum == 111
        bigger = find(trial_start_run2 > 183+84);
        for x = 1:length(bigger)
            trial_start_run2(bigger(x)) = trial_start_run2(bigger(x))+1;
        end
    end
        
        trials = horzcat(trial_start_run1, trial_start_run2);

% 
SF_hiconf = [];
for x = 1:length(trials)
    if mvpa_regs.preview.subresp(trials(x)) == 3
        SF_t = trials(x);
        SF_hiconf= horzcat(SF_hiconf, SF_t);
    elseif mvpa_regs.preview.subresp(trials(x)) == 4
        SF_t = trials(x);
        SF_hiconf= horzcat(SF_hiconf, SF_t);
    elseif mvpa_regs.preview.subresp(trials(x)) == 2
        SF_t = trials(x);
       SF_hiconf= horzcat(SF_hiconf, SF_t);
    else
        SF_hiconf = SF_hiconf;
    end
end

%
SR_hiconf = [];
for x = 1:length(trials)
    if mvpa_regs.preview.subresp(trials(x)) == 1
        SR_t = trials(x);
        SR_hiconf= horzcat(SR_hiconf, SR_t);
    else
        SR_hiconf = SR_hiconf;
    end
end

SF = [];
for x = 1:length(trials)
    if mvpa_regs.preview.subresp(trials(x)) == 3
        SF_t = trials(x);
        SF= horzcat(SF, SF_t);
    elseif mvpa_regs.preview.subresp(trials(x)) == 4
        SF_t = trials(x);
        SF= horzcat(SF, SF_t);
    else
        SF = SF;
    end
end

%
SR = [];
for x = 1:length(trials)
    if mvpa_regs.preview.subresp(trials(x)) == 1
        SR_t = trials(x);
        SR= horzcat(SR, SR_t);
    elseif mvpa_regs.preview.subresp(trials(x)) == 2
        SR_t = trials(x);
        SR= horzcat(SR, SR_t);
    else
        SR = SR;
    end
end

%%%%%%% MODELS %%%%%%%%%%
% Model 1 Main effects of memory instruction

% onsets{1} = forget;
% onsets{2} = remember;
% 
% % names = {'face', 'remember'};
% % durations = ones(length(onsets),1)';
% % durations= num2cell(durations);
% % 
% % outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_RFonsets_%s.mat', subjID, phase, subjID);
% % save(outfname,'names','onsets','durations');


% Model 2 SM effects of memory instruction

clear onsets
onsets{1,1} = SF_hiconf;
onsets{1,2} = SR_hiconf;


names = {'SF_hiconf', 'SR_hiconf',};
durations = ones(length(onsets),1)';
durations= num2cell(durations);

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_preview_RF_SM_hiconf2_onsets_%s.mat', subjID, subjID);
save(outfname,'names','onsets','durations');

clear onsets
onsets{1,1} = SF;
onsets{1,2} = SR;

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_preview_RF_SM2_onsets_%s.mat', subjID, subjID);

save(outfname,'names','onsets','durations');


%Model 3 emotional effects of memory instruction

% onsets{1} = for_neg;
% onsets{2} = for_neu;
% onsets{3} = rem_neg;
% onsets{4} = rem_neu;
% 
% names = {'for_neg', 'for_neu', 'rem_neg', 'rem_neu'};
% durations = ones(length(onsets),1)';
% durations= num2cell(durations);
% 
% outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_RF_emo_onsets_%s.mat', subjID, phase, subjID);
% save(outfname,'names','onsets','durations');
%

cd(script_dir);

end











