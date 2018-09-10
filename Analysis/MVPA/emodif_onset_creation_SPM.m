function emodif_onset_creation_SPM (subjNum)
%%
%this function converts the IMDiF onsets optimized for the Princeton MVPA toolbox
%for SPM - for both item onsets and instruction onsets.
%%
subjID = sprintf('emodif_%s',num2str(subjNum));
%subj_dir = sprintf('/corral-repl/utexas/lewpealab/imdif/%s', subjID);

%SET names for emodif_onsets to be read out

outfname_item = sprintf('~/emodif_data/emodif_%s/behav/EmoDif_SPM_item_onsets_%s.mat', subjID, subjID);
outfname_item = sprintf('~/emodif_data/emodif_%s/behav/EmoDif_SPM_instr_onsets_%s.mat', subjID, subjID);

%read in old onsets
load(sprintf('~/emodif_data/emodif_%s/behav/EmoDif_mvpa_allregs.mat', subjID));

%%%%% these variables ALL change with experiment %%%%
%write names and durations
localizer_run_TR = 342; %426 in subject 1
% preview_run_TR = 366;
% DFencode_run_TR = 426;

names = {'Face', 'Scene', 'Object', 'Word', 'Rest'};

% duration is 0 if we are modeling a delta - this goes up in number for
% boxcar models etc.  1 duration for each condition.

durations = {0 0 0 0 0};

%%%%%%%% BUILDING our onsets %%%%%%%%%%%%%%
% for each of these conditions:
%face scene object word rest
%%% Start Massive Indexing Effort or SMIE %%%

study_trials = EmoDif_mvpa_allregs.localizer.trial;
% trial num WITHIN a run (will expand this out later)

%%%%%% this is where i stopped. i can work on emodif_101 who has different
%%%%%% localizer regressors and we have to take a random sample between
%%%%%% neutral and negative - and the rest would go into noncrit regressor.
%%%%%% 

face_index = mvpa_regs.study.cat==1;
scene_index = mvpa_regs.study.cat==2;

% 1=Face 2=Scene

% run1_index = mvpa_regs.study.run == 1;
% run2_index = mvpa_regs.study.run == 2;
% run3_index = mvpa_regs.study.run == 3;
% run4_index = mvpa_regs.study.run == 4;
% run5_index = mvpa_regs.study.run == 5;
% run6_index = mvpa_regs.study.run == 6;
% % 1-6 Runs

rem_index = mvpa_regs.study.instr == 1;
for_index = mvpa_regs.study.instr == 2;

% 1= Remember 2 = Forget
c_old_index = mvpa_regs.test.resp == 1;
uc_old_index = mvpa_regs.test.resp == 2;
uc_new_index = mvpa_regs.test.resp == 3;
c_new_index = mvpa_regs.test.resp == 4;

old_index = c_old_index+uc_old_index;
new_index = c_new_index+uc_new_index;

%MAJOR FUCKUP WARNING 3-4 = subsequently forgotten, 1-2 subsequently remembered. 

%%% Start Defining Conditions by Variables from old onsets %%%%
%'fa_SF_FI'= a subsequently forgotten face trial with a forget instruction

fa_SF_FI_loc = find(face_index & new_index & for_index == 1);
fa_SR_FI_loc = find(face_index & old_index & for_index == 1);
fa_SF_RI_loc = find(face_index & new_index & rem_index == 1);
fa_SR_RI_loc = find(face_index & old_index & rem_index == 1);
sc_SF_FI_loc = find(scene_index & new_index & for_index == 1);
sc_SR_FI_loc = find(scene_index & old_index & for_index == 1);
sc_SF_RI_loc = find(scene_index & new_index & rem_index == 1);
sc_SR_RI_loc = find(scene_index & old_index & rem_index == 1);

%%% identified column numbers of these conditions %%%
%%% identify trial numbers and expand to find scan volume numbers that
%%% correspond with 1. the onset of the item and 2. the onset of the
%%% instruction.

% 1. fa_SF_FI

fa_SF_FI_items = [];
fa_SF_FI_instrs = [];

trial_index = study_trials(fa_SF_FI_loc);
run_index = mvpa_regs.study.run(fa_SF_FI_loc);

for i = 1:numel(trial_index)
    fa_SF_FI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    fa_SF_FI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    fa_SF_FI_items = horzcat(fa_SF_FI_items, fa_SF_FI_item);
    fa_SF_FI_instrs = horzcat(fa_SF_FI_instrs, fa_SF_FI_instr);
end

clear trial_index
clear run_index

% 2. fa_SR_FI

fa_SR_FI_items = [];
fa_SR_FI_instrs = [];

trial_index = study_trials(fa_SR_FI_loc);
run_index = mvpa_regs.study.run(fa_SR_FI_loc);

for i = 1:numel(trial_index)
    fa_SR_FI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    fa_SR_FI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    fa_SR_FI_items = horzcat(fa_SR_FI_items, fa_SR_FI_item);
    fa_SR_FI_instrs = horzcat(fa_SR_FI_instrs, fa_SR_FI_instr);
end

clear trial_index
clear run_index


% 3. fa_SF_RI

fa_SF_RI_items = [];
fa_SF_RI_instrs = [];

trial_index = study_trials(fa_SF_RI_loc);
run_index = mvpa_regs.study.run(fa_SF_RI_loc);

for i = 1:numel(trial_index)
    fa_SF_RI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    fa_SF_RI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    fa_SF_RI_items = horzcat(fa_SF_RI_items, fa_SF_RI_item);
    fa_SF_RI_instrs = horzcat(fa_SF_RI_instrs, fa_SF_RI_instr);
end

clear trial_index
clear run_index

% 4. fa_SR_RI

fa_SR_RI_items = [];
fa_SR_RI_instrs = [];

trial_index = study_trials(fa_SR_RI_loc);
run_index = mvpa_regs.study.run(fa_SR_RI_loc);

for i = 1:numel(trial_index)
    fa_SR_RI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    fa_SR_RI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    fa_SR_RI_items = horzcat(fa_SR_RI_items, fa_SR_RI_item);
    fa_SR_RI_instrs = horzcat(fa_SR_RI_instrs, fa_SR_RI_instr);
end

clear trial_index
clear run_index


% 5. sc_SF_FI

sc_SF_FI_items = [];
sc_SF_FI_instrs = [];

trial_index = study_trials(sc_SF_FI_loc);
run_index = mvpa_regs.study.run(sc_SF_FI_loc);

for i = 1:numel(trial_index)
    sc_SF_FI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    sc_SF_FI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    sc_SF_FI_items = horzcat(sc_SF_FI_items, sc_SF_FI_item);
    sc_SF_FI_instrs = horzcat(sc_SF_FI_instrs, sc_SF_FI_instr);
end

clear trial_index
clear run_index

% 6. fa_SR_FI

sc_SR_FI_items = [];
sc_SR_FI_instrs = [];

trial_index = study_trials(sc_SR_FI_loc);
run_index = mvpa_regs.study.run(sc_SR_FI_loc);

for i = 1:numel(trial_index)
    sc_SR_FI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    sc_SR_FI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    sc_SR_FI_items = horzcat(sc_SR_FI_items, sc_SR_FI_item);
    sc_SR_FI_instrs = horzcat(sc_SR_FI_instrs, sc_SR_FI_instr);
end

clear trial_index
clear run_index


% 7. sc_SF_RI

sc_SF_RI_items = [];
sc_SF_RI_instrs = [];

trial_index = study_trials(sc_SF_RI_loc);
run_index = mvpa_regs.study.run(sc_SF_RI_loc);

for i = 1:numel(trial_index)
    sc_SF_RI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    sc_SF_RI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    sc_SF_RI_items = horzcat(sc_SF_RI_items, sc_SF_RI_item);
    sc_SF_RI_instrs = horzcat(sc_SF_RI_instrs, sc_SF_RI_instr);
end

clear trial_index
clear run_index

% 8. sc_SR_RI

sc_SR_RI_items = [];
sc_SR_RI_instrs = [];

trial_index = study_trials(sc_SR_RI_loc);
run_index = mvpa_regs.study.run(sc_SR_RI_loc);

for i = 1:numel(trial_index)
    sc_SR_RI_item = ((trial_index(i)*9)-8)+(run_volnum*((run_index(i)-1)));
    sc_SR_RI_instr = ((trial_index(i)*9)-5)+(run_volnum*((run_index(i)-1)));
    sc_SR_RI_items = horzcat(sc_SR_RI_items, sc_SR_RI_item);
    sc_SR_RI_instrs = horzcat(sc_SR_RI_instrs, sc_SR_RI_instr);
end

clear trial_index
clear run_index

onsets{1} = fa_SF_FI_items;
onsets{2} = fa_SR_FI_items;
onsets{3} = fa_SF_RI_items;
onsets{4} = fa_SR_RI_items;
onsets{5} = sc_SF_FI_items;
onsets{6} = sc_SR_FI_items;
onsets{7} = sc_SF_RI_items;
onsets{8} = sc_SR_RI_items;

save(outfname_item,'names','onsets','durations');

clear onsets

onsets{1} = fa_SF_FI_instrs;
onsets{2} = fa_SR_FI_instrs;
onsets{3} = fa_SF_RI_instrs;
onsets{4} = fa_SR_RI_instrs;
onsets{5} = sc_SF_FI_instrs;
onsets{6} = sc_SR_FI_instrs;
onsets{7} = sc_SF_RI_instrs;
onsets{8} = sc_SR_RI_instrs;

save(outfname_instr,'names','onsets','durations');

clear all;

end
 
 