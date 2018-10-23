function emodif_onset_creation_SPM_hiconf (subjNum, phase)
%%
%this function converts the IMDiF onsets optimized for the Princeton MVPA toolbox
%for SPM - for all phases - Preview, DFencode_instr, DFencode_item,
%Localizer
%%
subjID = sprintf('emodif_%s',num2str(subjNum));
%subj_dir = sprintf('/corral-repl/utexas/lewpealab/imdif/%s', subjID);

%SET names for emodif_onsets to be read out

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_onsets_%s.mat', subjID, phase, subjID);
script_dir = pwd;

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



if strcmp(phase,'localizer') == 1
    
durations= repmat(9,length(names_local),1);
durations= num2cell(durations)';
    
if subjNum == 101 || subjNum == 102 || subjNum == 103
miniblock_run = 15;
else
miniblock_run= 12;
end


trial_start_run1 = 1:(loc_miniblock_TR+loc_postblock_rest):(loc_miniblock_TR+loc_postblock_rest)*miniblock_run;
trial_start_run2 = trial_start_run1(end)+ (loc_miniblock_TR+loc_postblock_rest) + betweenrun_rest:loc_miniblock_TR+loc_postblock_rest:(trial_start_run1(end)+betweenrun_rest)+((loc_miniblock_TR+loc_postblock_rest)*miniblock_run);


%here are exceptions for LOCALIZER
if subjNum == 101 %218 is 6 - extra TR
    bigger = find(trial_start_run2 > 218);
    for x = 1:length(bigger)
        trial_start_run2(bigger(x)) = trial_start_run2(bigger(x))+1;
    end
elseif subjNum == 113
        bigger = find(trial_start_run1 > 5);
    for x = 1:length(bigger)
        trial_start_run1(bigger(x)) = trial_start_run1(bigger(x))+1;
    end
elseif subjNum == 115
            bigger = find(trial_start_run2> (213+102));
    for x = 1:length(bigger)
        trial_start_run2(bigger(x)) = trial_start_run2(bigger(x))+1;
    end
elseif subjNum == 117
                bigger = find(trial_start_run1> 119);
    for x = 1:length(bigger)
        trial_start_run1(bigger(x)) = trial_start_run1(bigger(x))+1;
    end
elseif subjNum == 118 
                bigger = find(trial_start_run2> (213+132));
    for x = 1:length(bigger)
        trial_start_run2(bigger(x)) = trial_start_run2(bigger(x))+1;
    end
elseif subjNum == 124   
                bigger = find(trial_start_run1> 57);
    for x = 1:length(bigger)
        trial_start_run1(bigger(x)) = trial_start_run1(bigger(x))+1;
    end
end
trials = horzcat(trial_start_run1, trial_start_run2);


%for face blocks
face_blocks = [];

for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 1
        face_t = trials(x);
        face_blocks = horzcat(face_blocks,face_t);
    else
        face_blocks = face_blocks;
    end
end

%for scene blocks
scene_blocks = [];
for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 2
        scene_t = trials(x);
        scene_blocks = horzcat(scene_blocks,scene_t);
    else
        scene_blocks = scene_blocks;
    end
end
    
%for object blocks
object_blocks = [];
for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 3
        object_t = trials(x);
        object_blocks = horzcat(object_blocks,object_t);
    else
        object_blocks = object_blocks;
    end
end

%for word blocks
if subjNum == 101 || subjNum == 102 || subjNum == 103
    

neg_words_blocks = [];
for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 4
        neg_words_t = trials(x);
        neg_words_blocks = horzcat(neg_words_blocks,neg_words_t);
    else
        neg_words_blocks = neg_words_blocks;
    end
end

neu_words_blocks = [];
for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 5
        neu_words_t = trials(x);
        neu_words_blocks = horzcat(neu_words_blocks,neu_words_t);
    else
        neu_words_blocks = neu_words_blocks;
    end
end

neg_words_blocks = neg_words_blocks(randperm(length(neg_words_blocks)));
neu_words_blocks = neu_words_blocks(randperm(length(neu_words_blocks)));
word_blocks = sort(horzcat(neg_words_blocks(1:3),neu_words_blocks(1:3)));
noncritword_blocks = sort(horzcat(neg_words_blocks(4:6),neu_words_blocks(4:6)));

onsets{1} = face_blocks;
onsets{2} = scene_blocks;
onsets{3} = object_blocks;
onsets{4} = word_blocks;
onsets{5} = noncritword_blocks;

names = {'face', 'scene', 'object', 'word', 'noncritword'};

else
word_blocks = [];
for x = 1:length(trials)
    if mvpa_regs.localizer.cat(trials(x)) == 4
        words_t = trials(x);
       word_blocks = horzcat(word_blocks,words_t);
    else
        word_blocks = word_blocks;
    end   
end

onsets{1} = face_blocks;
onsets{2} = scene_blocks;
onsets{3} = object_blocks;
onsets{4} = word_blocks;
names = {'face', 'scene', 'object', 'word'};

save(outfname,'names','onsets','durations');

end

durations = durations(1:length(onsets));

elseif strcmp(phase,'DFEncode') == 1
    
    % three different models will be produced. 
    %model 1 Main effects of Memory instruction
    %model 2 SM effects - there will be some zero here - so use in between
    %rest for those conditions. 
    
    %model 3 Emotion effects
    
% names are going to be 
% SR_R subsequently remembered, remember instruction
% SR_F subsequently remembered, forget instruction
% SF_R subsequently forgotten, remember instruction
% SF_F subsequently forgotten, forget instruction

    trial_start_run1 = 1:DFEncode_trialTR:(DFEncode_trialTR*DFEncode_trialN)/2;
    trial_start_run2 = (DFEncode_trialTR*(DFEncode_trialN/2)+betweenrun_rest)+1:DFEncode_trialTR:(DFEncode_trialTR*DFEncode_trialN)+betweenrun_rest;
    
    %here are exceptions for DFencode
    
    if subjNum == 116 %
        bigger = find(trial_start_run2 > 213+162);
        for x = 1:length(bigger)
            trial_start_run2(bigger(x)) = trial_start_run2(bigger(x))+1;
        end
    end

trials = horzcat(trial_start_run1, trial_start_run2);


%% all vectors 
%instructions 0 = forget and 1 = remember
forget = [];
for x = 1:length(trials)
    if mvpa_regs.DFEncode.instr(trials(x)) == 0
        forget_t = trials(x);
        forget= horzcat(forget,forget_t);
    else
        forget =  forget;
    end
end

% for_neg

for_neg = [];
for x = 1:length(forget)
    if mvpa_regs.DFEncode.emo(forget(x)) == 1
        for_neg_t = forget(x);
        for_neg = horzcat(for_neg, for_neg_t);
    else
        for_neg = for_neg;
    end
end

% 
for_neu = [];
for x = 1:length(forget)
    if mvpa_regs.DFEncode.emo(forget(x)) == 0
        for_neu_t = forget(x);
        for_neu = horzcat(for_neu, for_neu_t);
    else
        for_neu = for_neu;
    end
end 

% 
for_SF = [];
for x = 1:length(forget)
    if mvpa_regs.DFEncode.subresp(forget(x)) == 3
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    elseif mvpa_regs.DFEncode.subresp(forget(x)) == 4
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    elseif mvpa_regs.DFEncode.subresp(forget(x)) == 2
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    else
        for_SF = for_SF;
    end
end

% 
for_SR = [];
for x = 1:length(forget)
    if mvpa_regs.DFEncode.subresp(forget(x)) == 1
        for_SR_t = forget(x);
        for_SR= horzcat(for_SR, for_SR_t);
    else
        for_SR = for_SR;
    end
end

remember = [];
for x = 1:length(trials)
    if mvpa_regs.DFEncode.instr(trials(x)) == 1
        remember_t = trials(x);
        remember= horzcat(remember,remember_t);
    else
        remember =  remember;
    end
end

% for_neg

rem_neg = [];
for x = 1:length(remember)
    if mvpa_regs.DFEncode.emo(remember(x)) == 1
        rem_neg_t = remember(x);
        rem_neg = horzcat(rem_neg, rem_neg_t);
    else
        rem_neg = rem_neg;
    end
end

% 
rem_neu = [];
for x = 1:length(remember)
    if mvpa_regs.DFEncode.emo(remember(x)) == 0
        rem_neu_t = remember(x);
        rem_neu = horzcat(rem_neu, rem_neu_t);
    else
        rem_neu = rem_neu;
    end
end 

% 
rem_SF = [];
for x = 1:length(remember)
    if mvpa_regs.DFEncode.subresp(remember(x)) == 3
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    elseif mvpa_regs.DFEncode.subresp(remember(x)) == 4
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    elseif mvpa_regs.DFEncode.subresp(remember(x)) == 2
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    else
        rem_SF = rem_SF;
    end
end

% 
rem_SR = [];
for x = 1:length(remember)
    if mvpa_regs.DFEncode.subresp(remember(x)) == 1
        rem_SR_t = remember(x);
        rem_SR= horzcat(rem_SR, rem_SR_t);
    else
        rem_SR = rem_SR;
    end
end

%%%%%%% MODELS %%%%%%%%%%
% Model 1 Main effects of memory instruction
% 
% onsets{1} = forget;
% onsets{2} = remember;
% 
% names = {'forget', 'remember'};
% durations = ones(length(onsets),1)';
% durations= num2cell(durations);
% 
% outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_RFonsets_%s.mat', subjID, phase, subjID);
% save(outfname,'names','onsets','durations');

% Model 2 SM effects of memory instruction

onsets{1} = for_SF;
onsets{2} = for_SR;
onsets{3} = rem_SF;
onsets{4} = rem_SR;

names = {'for_SF', 'for_SR', 'rem_SF', 'rem_SR'};
durations = ones(length(onsets),1)';
durations= num2cell(durations);

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_RF_SM_hiconf_onsets_%s.mat', subjID, phase, subjID);
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

elseif strcmp(phase,'preview') ==  1
    
    % three different models will be produced. 
    %model 1 Main effects of Memory instruction
    %model 2 SM effects - there will be some zero here - so use in between
    %rest for those conditions. 
    
    %model 3 Emotion effects
    
    %
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
        
        %% all vectors 
%instructions 0 = forget and 1 = remember
forget = [];
for x = 1:length(trials)
    if mvpa_regs.preview.instr(trials(x)) == 0
        forget_t = trials(x);
        forget= horzcat(forget,forget_t);
    else
        forget =  forget;
    end
end

% for_neg

for_neg = [];
for x = 1:length(forget)
    if mvpa_regs.preview.emo(forget(x)) == 1
        for_neg_t = forget(x);
        for_neg = horzcat(for_neg, for_neg_t);
    else
        for_neg = for_neg;
    end
end

% 
for_neu = [];
for x = 1:length(forget)
    if mvpa_regs.preview.emo(forget(x)) == 0
        for_neu_t = forget(x);
        for_neu = horzcat(for_neu, for_neu_t);
    else
        for_neu = for_neu;
    end
end 

% 
for_SF = [];
for x = 1:length(forget)
    if mvpa_regs.preview.subresp(forget(x)) == 3
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    elseif mvpa_regs.preview.subresp(forget(x)) == 4
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    elseif mvpa_regs.preview.subresp(forget(x)) == 2
        for_SF_t = forget(x);
        for_SF= horzcat(for_SF, for_SF_t);
    else
        for_SF = for_SF;
    end
end

% 
for_SR = [];
for x = 1:length(forget)
    if mvpa_regs.preview.subresp(forget(x)) == 1
        for_SR_t = forget(x);
        for_SR= horzcat(for_SR, for_SR_t);
    else
        for_SR = for_SR;
    end
end

remember = [];
for x = 1:length(trials)
    if mvpa_regs.preview.instr(trials(x)) == 1
        remember_t = trials(x);
        remember= horzcat(remember,remember_t);
    else
        remember =  remember;
    end
end

% for_neg

rem_neg = [];
for x = 1:length(remember)
    if mvpa_regs.preview.emo(remember(x)) == 1
        rem_neg_t = remember(x);
        rem_neg = horzcat(rem_neg, rem_neg_t);
    else
        rem_neg = rem_neg;
    end
end

% 
rem_neu = [];
for x = 1:length(remember)
    if mvpa_regs.preview.emo(remember(x)) == 0
        rem_neu_t = remember(x);
        rem_neu = horzcat(rem_neu, rem_neu_t);
    else
        rem_neu = rem_neu;
    end
end 

% 
rem_SF = [];
for x = 1:length(remember)
    if mvpa_regs.preview.subresp(remember(x)) == 3
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    elseif mvpa_regs.preview.subresp(remember(x)) == 4
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    elseif mvpa_regs.preview.subresp(remember(x)) == 2
        rem_SF_t = remember(x);
        rem_SF= horzcat(rem_SF, rem_SF_t);
    else
        rem_SF = rem_SF;
    end
end

% 
rem_SR = [];
for x = 1:length(remember)
    if mvpa_regs.preview.subresp(remember(x)) == 1
        rem_SR_t = remember(x);
        rem_SR= horzcat(rem_SR, rem_SR_t);
    else
        rem_SR = rem_SR;
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

onsets{1} = for_SF;
onsets{2} = for_SR;
onsets{3} = rem_SF;
onsets{4} = rem_SR;

names = {'for_SF', 'for_SR', 'rem_SF', 'rem_SR'};
durations = ones(length(onsets),1)';
durations= num2cell(durations);

outfname = sprintf('~/emodif_data/%s/behav/EmoDif_SPM_%s_RF_SM_hiconf_onsets_%s.mat', subjID, phase, subjID);
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

end

cd(script_dir);

end
 
 