function [rsa] = emodif_rsa_aggL_P_DF_bycue(subjNum,maskName, shift, dfencodeTRstart, dfencodeTRlength)
%  [rsa] = emodif_rsa_aggL_P_DF_bycue('101','Scene100', 2, 2, 3)
%RSA set up script based off of hyojeong's clearmem matlab RSA script. ***
%requires princeton toolbox ***

%this RSA is to compare PREVIEW scene related activity with ENCODING scene
%related activity following the same word presentation. 

%%%% THIS VERSION AGGREGATES RSA Localizer BEFORE correlation %%%%

%shift is TR shift for training.

%Experiment:EmoDF
%Type : Emotional Directed Forgetting with Scene tags
%Phase: Preview to Study
% Date      : August 2018
% Version   : 1
% Author    : Tracy Wang
% Contact   : tracy.wang@utexas.edu


  
  %%% parameters %%%
  
  start_dir = pwd;
  args.experiment = 'emodif';
  args.train_phase = 'preview';
  args.test_phase = 'DFencode';
  args.preview.nTRs = 366;
  args.DFencode.nTRs = 426;
  if subjNum == '101' 
      args.Localizer.nTRs = 426;
  elseif subjNum == '102'
      args.Localizer.nTRs = 426;
  elseif subjNum == '103'
      args.Localizer.nTRs = 426;
  else
      args.Localizer.nTRs = 342; %426 in subject 1, 2 and 3.  %342 for everyone else.
  end
 
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  args.preview.trial_length = 6;
  args.DFencode.trial_length = 7;
  args.Localizer.trial_length = 14; %9 TRs of miniblocks, 5 TRs of rest
  args.preview.trial_break = 3;
  args.DFencode.trial_break = 3;
  args.Localizer.trial_break = 3;
  args.preview.trialnum = 60;
  args.DFencode.trialnum = 60;
    if subjNum == '101'
      args.Localizer.trialnum = 30;
  elseif subjNum == '102'
      args.Localizer.trialnum = 30;
  elseif subjNum == '103'
      args.Localizer.trialnum = 30;
      
  else
      args.Localizer.trialnum = 24; %426 in subject 1, 2 and 3.  %342 for everyone else.
  end
  %adjustable parameters
%%%% now doing both %%%%
  %1 = mean (preview) x mean (DFencode) uses meanTR_length
  %2 = mean (preview) x TR (DF encode) uses only meanTR_length of preview
  
  args.preview.meanTR_length = 3;
  args.preview.meanTR_start = 3;
  args.Localizer.meanTR_start = 3;
  args.Localizer.meanTR_length = 6;
  args.DFencode.meanTR_length = dfencodeTRlength; % last TR goes into DF instruction, so could be 3TRs, must try empirically
  args.DFencode.meanTR_start = dfencodeTRstart;
  args.preview.meanTR_end = (args.preview.meanTR_start+args.preview.meanTR_length)-1;
  args.DFencode.meanTR_end = (args.DFencode.meanTR_start+args.DFencode.meanTR_length)-1;
  args.Localizer.meanTR_end = (args.Localizer.meanTR_start+args.Localizer.meanTR_length)-1;
  
  args.shiftTR = shift; %TR train shift. 
  
  %for astoria
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
  %for tigger
  %args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
  
  %data directories
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
%   args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.mask_dir = sprintf('%s/roi_generation/SPM_GLM', args.subj_dir)
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
%   args.DFencode_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.test_phase, test_date);
%   args.preview_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.train_phase, train_date);
  args.output_dir = sprintf('%s/results/rsa_results/preview_DF_preview/%s/',args.subj_dir,maskName);
  mkdir(args.output_dir);
  args.subjNum = subjNum;
  
  %----------------------------------------------------------------------
  % turn on diary to capture analysis output
  %
  diary on;
  diary(sprintf('%s/%s_diary.txt',args.output_dir,args.subjNum));
  fprintf('###########################################\n\n');
  disp(args);
  fprintf('###########################################\n');
  %
%%%% LOAD REGRESSORS & RSA Parameters %

load(sprintf('%s/EmoDiF_mvpa_allregs.mat', args.regs_dir));
load('/Users/tw24955/emodif_data/emodif_counterbalance.mat');


if subjNum == '101' | subjNum =='102'| subjNum == '103'| subjNum == '104'

load(sprintf('/Users/tw24955/emodif_data/RSA_params.mat')); % just the same for the first 4 subjects 

%%%% 

%%% MANUALLY adjust regressors <note, this is where we would spike correct
%%% by zeroing out those TRs with bad spikes%%% currently nothing here

%%% expanding RSA parameters with mvpa_regs %%%


rsa.preview.preview2DFencode = RSA_params.preview2study;
rsa.preview.preview2DFencode_nonexpanded = RSA_params.NonExpandedPrv2Sty;
rsa.DFencode.DFencode2preview_nonexpanded = RSA_params.NonExpandedSty2Prv;
rsa.DFencode.DFencode2preview = RSA_params.study2preview;
else
rsa.preview.preview2DFencode = mvpa_regs.preview.ExpandedStudyTrialMapping;
rsa.preview.preview2DFencode_nonexpanded = mvpa_regs.preview.StudyTrialMapping;
rsa.DFencode.DFencode2preview_nonexpanded = mvpa_regs.DFEncode.PreviewTrialMapping;  
rsa.DFencode.DFencode2preview = mvpa_regs.DFEncode.ExpandedPreviewTrialMapping;
    
end

rsa.DFencode.instr = mvpa_regs.DFEncode.instr;

%preview trial numbers

firstblocktrialnum = [];
secondblocktrialnum = [];
for i = 1:(args.preview.trialnum/2)
    x = repmat(i,args.preview.trial_length,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for j = ((args.preview.trialnum/2)+1):args.preview.trialnum
    x = repmat(j,args.preview.trial_length,1)';
    secondblocktrialnum = horzcat(secondblocktrialnum, x);
end

previewtrialbreak = zeros(args.preview.trial_break,1)';

rsa.preview.trialnum = horzcat(firstblocktrialnum,previewtrialbreak,secondblocktrialnum,previewtrialbreak);

%DFencode trial numbers

firstblocktrialnum = [];
secondblocktrialnum = [];

for i = 1:30
    x = repmat(i,7,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for i = 31:60
    x = repmat(i,7,1)';
    secondblocktrialnum = horzcat(secondblocktrialnum, x);
end

DFencodetrialbreak = zeros(args.DFencode.trial_break,1)';

rsa.DFencode.trialnum = horzcat(firstblocktrialnum,previewtrialbreak,secondblocktrialnum,DFencodetrialbreak);

%Localizer trial numbers

firstblocktrialnum = [];
secondblocktrialnum = [];

for i = 1:(args.Localizer.trialnum/2)
    x = repmat(i,args.Localizer.trial_length,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for i = ((args.Localizer.trialnum/2)+1):args.Localizer.trialnum
    x = repmat(i,args.Localizer.trial_length,1)';
    secondblocktrialnum = horzcat(secondblocktrialnum, x);
end

Localizertrialbreak = zeros(args.Localizer.trial_break,1)';

rsa.Localizer.trialnum = horzcat(firstblocktrialnum,Localizertrialbreak,secondblocktrialnum,Localizertrialbreak);

if subjNum == '101'
    trialdummy = horzcat(0,rsa.Localizer.trialnum);
    rsa.Localizer.trialnum(218:end) = trialdummy(218:length(rsa.Localizer.trialnum));

end

% bring DFencode category and responses into preview regressors
for i = 1:args.preview.nTRs
    if rsa.preview.preview2DFencode(i) == 0
        trial_cat = NaN;
    else trial_dfencode_num = find(mvpa_regs.DFEncode.trial == rsa.preview.preview2DFencode(i));
     trial_cat = mvpa_regs.DFEncode.cat(1,trial_dfencode_num(1));
     trial_resp = mvpa_regs.DFEncode.subresp(1,trial_dfencode_num(1));
    end
    rsa.preview.DFencodecat(i) = trial_cat;
    rsa.preview.DFencodeResp(i) = trial_resp;
end

%complete DFencode regressors for RSA


for i = 1:args.DFencode.nTRs
    if rsa.DFencode.DFencode2preview(i) == 0
        trial_cat = NaN;
    else trial_DFencode = find(mvpa_regs.preview.trial == rsa.DFencode.DFencode2preview(i));
        trial_cat = mvpa_regs.DFEncode.cat(1,trial_DFencode(1));
        trial_resp = mvpa_regs.DFEncode.subresp(1,trial_DFencode(1));
        
    end
    rsa.DFencode.cat(i) = trial_cat;
    rsa.DFencode.resp(i) = trial_resp;
end
        
% directories for across subject RSA
subjNum_dbl = str2num(subjNum);
subjNum_position = find(counterbalance(:,1)==subjNum_dbl);
subjNum2 = num2str(counterbalance(subjNum_position,2));
subjNum3 = num2str(counterbalance(subjNum_position,3));
subjNum4 = num2str(counterbalance(subjNum_position,4));

args.subjID2 = sprintf('emodif_%s',num2str(subjNum2));
args.subjID3 = sprintf('emodif_%s',num2str(subjNum3));
args.subjID4 = sprintf('emodif_%s',num2str(subjNum4));

args.subj_dir2 = sprintf('/Users/tw24955/emodif_data/%s', args.subjID2);
args.subj_dir3 = sprintf('/Users/tw24955/emodif_data/%s', args.subjID3);
args.subj_dir4= sprintf('/Users/tw24955/emodif_data/%s', args.subjID4);

args.bold_dir2 = sprintf('%s/BOLD', args.subj_dir2);
args.bold_dir3 = sprintf('%s/BOLD', args.subj_dir3);
args.bold_dir4 = sprintf('%s/BOLD', args.subj_dir4);

args.mask_dir2 = sprintf('%s/roi_generation/SPM_GLM', args.subj_dir2);
args.mask_dir3 = sprintf('%s/roi_generation/SPM_GLM', args.subj_dir3);
args.mask_dir4 = sprintf('%s/roi_generation/SPM_GLM', args.subj_dir4);

%load masks - also need masks from 3 reciprocal subjects. 

mvpa_mask = fullfile(args.mask_dir, sprintf('%s.nii',maskName));
mvpa_mask2 = fullfile(args.mask_dir2, sprintf('%s.nii',maskName));
mvpa_mask3 = fullfile(args.mask_dir3, sprintf('%s.nii',maskName));
mvpa_mask4 = fullfile(args.mask_dir4, sprintf('%s.nii',maskName));

%% ============= Initializing subj. structure:start by creating an empty subj structure
% summarize(subj): summarize all info in subj structure
% get_object/set_object/set_objfield/set_objsubfield
% get_mat/set_mat

subj = init_subj(args.experiment, args.subjID);%identifier of the subj
fprintf('\n(+) %s phase data\n\n', args.train_phase);


%% ============= 01: EPI PATTERNS
%*************** load mask + read in epis

subj = load_spm_mask(subj,maskName,mvpa_mask);

mask_voxel = get_mat(subj,'mask', maskName);

fprintf('\n(+) load epi data under mask with %s voxels\n', num2str(count(mask_voxel)));

subj = load_spm_mask(subj,'mask_counter2',mvpa_mask2);
mask_voxel2 = get_mat(subj,'mask', 'mask_counter2');
fprintf('\n(+) load epi data under mask with %s voxels\n', num2str(count(mask_voxel2)));

subj = load_spm_mask(subj,'mask_counter3',mvpa_mask3);
mask_voxel3 = get_mat(subj,'mask', 'mask_counter3');
fprintf('\n(+) load epi data under mask with %s voxels\n', num2str(count(mask_voxel3)));

subj = load_spm_mask(subj,'mask_counter4',mvpa_mask4);
mask_voxel4 = get_mat(subj,'mask', 'mask_counter4');
fprintf('\n(+) load epi data under mask with %s voxels\n', num2str(count(mask_voxel4)));

Preview_raw_filenames = {'Preview1_corr_dt_mcf_brain.nii','Preview2_corr_dt_mcf_brain.nii'};
DFencode_raw_filenames={'DF_encoding_1_corr_dt_mcf_brain.nii', 'DF_encoding_2_corr_dt_mcf_brain.nii'};
Localizer_raw_filenames = {'MVPA_training_1_corr_dt_mcf_brain.nii', 'MVPA_training_2_corr_dt_mcf_brain.nii'};

cd(args.bold_dir)

subj = load_spm_pattern(subj, 'Preview_epis', maskName, Preview_raw_filenames);
subj = load_spm_pattern(subj, 'DFencode_epis', maskName, DFencode_raw_filenames);
subj = load_spm_pattern(subj, 'Localizer_epis', maskName, Localizer_raw_filenames);

cd(args.bold_dir2)
subj = load_spm_pattern(subj, 'Counter2_epis', 'mask_counter2', Preview_raw_filenames);

cd(args.bold_dir3)
subj = load_spm_pattern(subj, 'Counter3_epis', 'mask_counter3', Preview_raw_filenames);

cd(args.bold_dir4)
subj = load_spm_pattern(subj, 'Counter4_epis', 'mask_counter4', Preview_raw_filenames);

summarize(subj)

%*************** zsoring epis

%create selector to zscore across ALL runs - use full 'runs' selector

for x = 1:length(rsa.preview.trialnum);
    if x > args.preview.nTRs/2 
        all_preview_runs_selector(x) = 2;
    else all_preview_runs_selector(x) = 1;
    end
end

subj = initset_object(subj, 'selector', 'all_preview_runs', all_preview_runs_selector);

fprintf('[+] z-scoring Preview data\n\n');
subj = zscore_runs(subj,'Preview_epis', 'all_preview_runs');

%do the same for the other Counter runs
subj = zscore_runs(subj,'Counter2_epis', 'all_preview_runs');
subj = zscore_runs(subj,'Counter3_epis', 'all_preview_runs');
subj = zscore_runs(subj,'Counter4_epis', 'all_preview_runs');

for x = 1:length(rsa.DFencode.trialnum);
    if x >   args.DFencode.nTRs/2
        all_DFencode_runs_selector(x) = 2;
    else all_DFencode_runs_selector(x) = 1;
    end
end

subj = initset_object(subj, 'selector', 'all_DFencode_runs', all_DFencode_runs_selector);

  fprintf('[+] z-scoring DFencode data\n\n');
  subj = zscore_runs(subj,'DFencode_epis', 'all_DFencode_runs');
  
  for x = 1:length(rsa.Localizer.trialnum);
    if x >   args.Localizer.nTRs/2
        all_Localizer_runs_selector(x) = 2;
    else all_Localizer_runs_selector(x) = 1;
    end
end


subj = initset_object(subj, 'selector', 'all_Localizer_runs', all_Localizer_runs_selector);

fprintf('[+] z-scoring Localizer data\n\n');
subj = zscore_runs(subj,'Localizer_epis', 'all_Localizer_runs');



%*************** option: read epis in mask | wholebrain - NOT DONE
%currently everything is done through mask
 %*************** define: args.wholebrain '0 | 1' - NOT DONE
 
 
%% ============= Shift regressors and average pattern


%% ============== Average Patterns for option 1 (rsa_type)
%  if rsa_type = 1
%      
%      % take regressors and use training data to shift 2-3 TRs over
%  else

%averaging Localizer phase to ONLY scenes - 6 trials and the correlating
%separatly to Preview and DFencoding


for m = 1:args.Localizer.trialnum
    Localizer_trialTR_idx = find(rsa.Localizer.trialnum == m);
    Localizer_trialCAT = (mvpa_regs.localizer.cat(Localizer_trialTR_idx(dfencodeTRstart)));
    Localizer_trialTR_mean = Localizer_trialTR_idx(args.Localizer.meanTR_start:(args.Localizer.meanTR_start+args.Localizer.meanTR_length-1));
    Localizer_trial_patt = mean(subj.patterns{1,12}.mat(:,Localizer_trialTR_mean(1):(Localizer_trialTR_mean(args.Localizer.meanTR_length))),2);
    rsa.Localizer.mean.patterns(:,m) = Localizer_trial_patt;
    rsa.Localizer.mean.category(:,m) = Localizer_trialCAT;
end

    rsa.Localizer.mean.patterns_tosort = vertcat(rsa.Localizer.mean.category,rsa.Localizer.mean.patterns)';
    rsa.Localizer.mean.patterns_sorted = sortrows(rsa.Localizer.mean.patterns_tosort)';
    rsa.Localizer.mean.patterns_bycategory = rsa.Localizer.mean.patterns_sorted(2:end, :);
    rsa.Localizer.mean.patterns_scenes = rsa.Localizer.mean.patterns_sorted(2:end, 7:12);
    
    %Preview Phase
       
    for m = 1:args.preview.trialnum  
        trial_preview_idx = find(rsa.preview.trialnum == m);
        trialTR_mean = trial_preview_idx(args.preview.meanTR_start:(args.preview.meanTR_start+args.preview.meanTR_length-1));
        trial_patt = mean(subj.patterns{1,7}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        rsa.preview.orig.patterns(:,m) = trial_patt;
        
        counter2_patt = mean(subj.patterns{1,8}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        rsa.preview.counter2.patterns(:,m) = counter2_patt;
        
        counter3_patt = mean(subj.patterns{1,9}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        rsa.preview.counter3.patterns(:,m) = counter3_patt;
        
        counter4_patt = mean(subj.patterns{1,10}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        rsa.preview.counter4.patterns(:,m) = counter4_patt;
        
        rsa.preview.mean.patterns{m} = horzcat(rsa.preview.orig.patterns(:,m),rsa.preview.counter2.patterns(:,m), rsa.preview.counter3.patterns(:,m), rsa.preview.counter4.patterns(:,m));
        rsa.preview.mean.meanpattern(:,m) = mean(rsa.preview.mean.patterns{m},2); 
   
    end
    
    %matching to DFpatterns
    for m = 1:args.DFencode.trialnum
        %         trialTR_idx = find(rsa.DFencode.trialnum == m);
        %         %for shift
        %         trialTR_idx_sh = trialTR_idx+args.previewshiftTR;
        %         trialTR_mean = trialTR_idx_sh(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
        %         trial_patt = mean(subj.patterns{1,4}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        %         rsa.DFencode.mean.patterns(:,m) = trial_patt;
        
        P_trial_idx_match = find(rsa.preview.preview2DFencode_nonexpanded == m);
        P_trial_patt_mean_match = rsa.preview.mean.meanpattern(:,P_trial_idx_match);
        rsa.preview.mean.patterns_match(:,P_trial_idx_match) = P_trial_patt_mean_match;
    end

  %DFencode Phase
  
      for m = 1:args.DFencode.trialnum  
        trial_DFencode_idx = find(rsa.DFencode.trialnum == m);
        trialTR_mean = trial_DFencode_idx(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
        trial_patt = mean(subj.patterns{1,11}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
        rsa.DFencode.mean.patterns(:,m) = trial_patt;
      end
      
      
    
      % R and F
rsa.DFencode.mean.Rpatterns = [];
rsa.DFencode.mean.Fpatterns = [];
rsa.preview.mean.Fpatterns = [];
rsa.preview.mean.Rpatterns = [];


for k = 1:args.DFencode.trialnum
    DF_trialTR_idx = find(rsa.DFencode.trialnum == k);
    if rsa.DFencode.instr(DF_trialTR_idx(3)) == 0 %if the 1st TR shows that it's Forget
        %for shift
        DF_trialTR_Fidx_TRsOI = DF_trialTR_idx(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1)); %TRs of interest
        DF_trial_patt_Fmean = mean(subj.patterns{1,11}.mat(:,(DF_trialTR_Fidx_TRsOI(1)):DF_trialTR_Fidx_TRsOI(args.DFencode.meanTR_length)),2);
        rsa.DFencode.mean.Fpatterns = horzcat(rsa.DFencode.mean.Fpatterns,DF_trial_patt_Fmean);  
        %for preview
        P_trial_Fidx_match = find(rsa.preview.preview2DFencode_nonexpanded == k);
        P_trial_Fpatt_mean_match = rsa.preview.mean.meanpattern(:,P_trial_Fidx_match);
        rsa.preview.mean.Fpatterns = horzcat(rsa.preview.mean.Fpatterns,P_trial_Fpatt_mean_match);
        
    elseif rsa.DFencode.instr(DF_trialTR_idx(3)) == 1
        DF_trialTR_Ridx_TRsOI = DF_trialTR_idx(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1)); %TRs of interest
        DF_trial_patt_Rmean = mean(subj.patterns{1,11}.mat(:,(DF_trialTR_Ridx_TRsOI(1)):DF_trialTR_Ridx_TRsOI(args.DFencode.meanTR_length)),2);
        rsa.DFencode.mean.Rpatterns = horzcat(rsa.DFencode.mean.Rpatterns,DF_trial_patt_Rmean);
        
                %for preview
        P_trial_Ridx_match = find(rsa.preview.preview2DFencode_nonexpanded == k);
        P_trial_Rpatt_mean_match = rsa.preview.mean.meanpattern(:,P_trial_Ridx_match);
        rsa.preview.mean.Rpatterns = horzcat(rsa.preview.mean.Rpatterns,P_trial_Rpatt_mean_match);        
        
    end
end

%correlations


    corr_matrix_match_full = zeros(args.DFencode.trialnum,args.DFencode.trialnum);
    corr_matrix_match_fullz= zeros(args.DFencode.trialnum,args.DFencode.trialnum);
    
     for x = 1:args.DFencode.trialnum %trial number
        
        for y = 1:args.DFencode.trialnum %trial number
            
            corr_matrix_match_full(x, y) = corr2(rsa.DFencode.mean.patterns(:,x), rsa.preview.mean.patterns_match(:,y));
            corr_matrix_match_fullz(x, y) = 0.5*log((1+corr_matrix_match_full(x, y))/(1-corr_matrix_match_full(x, y)));
            rsa.results.smatrix.corr_matrix_match_full = corr_matrix_match_full;
            rsa.results.smatrix.corr_matrix_match_fullz = corr_matrix_match_fullz;
        end
     end
    
      corr_matrix_match_F = zeros(args.DFencode.trialnum/2,args.DFencode.trialnum/2);
    corr_matrix_match_Fz = zeros(args.DFencode.trialnum/2,args.DFencode.trialnum/2);
    
    for x = 1:args.DFencode.trialnum/2 %trial number
        
        for y = 1:args.DFencode.trialnum/2 %trial number
            
            corr_matrix_match_F(x, y) = corr2(rsa.DFencode.mean.Fpatterns(:,x), rsa.preview.mean.Fpatterns(:,y));
            corr_matrix_match_Fz(x, y) = 0.5*log((1+corr_matrix_match_F(x, y))/(1-corr_matrix_match_F(x, y)));
            
            rsa.results.smatrix.corr_matrix_match_F = corr_matrix_match_F;
            rsa.results.smatrix.corr_matrix_match_Fz = corr_matrix_match_Fz;
        end
    end
    
    corr_matrix_match_R = zeros(args.DFencode.trialnum/2,args.DFencode.trialnum/2);
    corr_matrix_match_Rz = zeros(args.DFencode.trialnum/2,args.DFencode.trialnum/2);
    
    for x = 1:args.DFencode.trialnum/2 %trial number
        
        for y = 1:args.DFencode.trialnum/2 %trial number
            
            corr_matrix_match_R(x, y) = corr2(rsa.DFencode.mean.Rpatterns(:,x), rsa.preview.mean.Rpatterns(:,y));
            corr_matrix_match_Rz(x, y) = 0.5*log((1+corr_matrix_match_R(x, y))/(1-corr_matrix_match_R(x, y)));
            rsa.results.smatrix.corr_matrix_match_R = corr_matrix_match_R;
            rsa.results.smatrix.corr_matrix_match_Rz = corr_matrix_match_Rz;
        end
    end
    
    cd(args.output_dir)
    
    fig = figure;
    set(fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    
    imagesc(corr_matrix_match_full); colormap('parula'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Preview agg patterns averaged over TR %d: TR %d run',args.DFencode.meanTR_start,args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode raw agg patterns averaged over TR %d: TR %d run',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(fig, sprintf('emodif%s_run_agg_full_TRstart%d_mean%d',subjNum,args.DFencode.meanTR_start,args.DFencode.meanTR_length),'png')    
    
        F_fig = figure;
    set(F_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    
    imagesc(corr_matrix_match_F); colormap('parula'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Preview agg patterns averaged over TR %d: TR %d run',args.DFencode.meanTR_start,args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode raw F patterns averaged over TR %d: TR %d run',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(F_fig, sprintf('emodif%s_run_agg_F_TRstart%d_mean%d',subjNum,args.DFencode.meanTR_start,args.DFencode.meanTR_length),'png') 
    
    
        R_fig = figure;
    set(R_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    
    imagesc(corr_matrix_match_R); colormap('parula'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Preview agg patterns averaged over TR %d: TR %d run',args.DFencode.meanTR_start,args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode raw R patterns averaged over TR %d: TR %d run',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(R_fig, sprintf('emodif%s_run_agg_R_TRstart%d_mean%d',subjNum,args.DFencode.meanTR_start,args.DFencode.meanTR_length),'png') 
    
    
    
    rsa.parameters = args;
    results_file = sprintf('%s/%s_TR%dto%d_Pagg_rsa_results.mat',...
           args.output_dir, args.subjID, args.DFencode.meanTR_start, args.DFencode.meanTR_end);
       save(results_file,'rsa')
       
    cd (start_dir); 
end
    
    
    
    
    
