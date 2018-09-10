function [rsa] = emodif_rsa_allphases(subjNum,maskName, train_date, test_date, shift)
%  [rsa] = emodif_rsa_preview_dfencode('103','tempoccfusi_pHg_LOC_combined_epi_space','29-Aug-2018', '29-Aug-2018', 2)
%RSA set up script based off of hyojeong's clearmem matlab RSA script. ***
%requires princeton toolbox ***

%this RSA is to compare PREVIEW scene related activity with ENCODING scene
%related activity following the same word presentation. 

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
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  args.preview.trial_length = 6;
  args.DFencode.trial_length = 7;
  args.preview.trial_break = 3;
  args.DFencode.trial_break = 3;
  args.trialnum = 60;
  
  %adjustable parameters
%%%% now doing both %%%%
  %1 = mean (preview) x mean (DFencode) uses meanTR_length
  %2 = mean (preview) x TR (DF encode) uses only meanTR_length of preview
  
  args.preview.meanTR_length = 3;
  args.preview.meanTR_start = 3;
  args.DFencode.meanTR_length = 3; % last TR goes into DF instruction, so could be 3TRs, must try empirically
  args.DFencode.meanTR_start = 2;
  args.shiftTR = shift; %TR train shift. 
  args.preview.meanTR_end = (args.preview.meanTR_start+args.preview.meanTR_length)-1;
  args.DFencode.meanTR_end = (args.DFencode.meanTR_start+args.DFencode.meanTR_length)-1;
  
  %for astoria
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
  %for tigger
  %args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
  
  %data directories
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.DFencode_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.test_phase, test_date);
  args.preview_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.train_phase, train_date);
  args.output_dir = sprintf('%s/results/rsa_results/%s/%s',args.subj_dir,maskName,shift);
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
load(sprintf('/Users/tw24955/emodif_data/RSA_params.mat')); % just the same for the first 4 subjects

%%%% 

%%% MANUALLY adjust regressors <note, this is where we would spike correct
%%% by zeroing out those TRs with bad spikes%%% currently nothing here

%%% expanding RSA parameters with mvpa_regs %%%

rsa.preview.preview2DFencode = RSA_params.preview2study;
rsa.preview.preview2DFencode_nonexpanded = RSA_params.NonExpandedPrv2Sty;
rsa.DFencode.DFencode2preview_nonexpanded = RSA_params.NonExpandedSty2Prv;

%preview trial numbers

firstblocktrialnum = [];
secondblocktrialnum = [];
for i = 1:(args.trialnum/2)
    x = repmat(i,args.preview.trial_length,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for j = ((args.trialnum/2)+1):args.trialnum
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
rsa.DFencode.DFencode2preview = RSA_params.study2preview;

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
        
%load masks

mvpa_mask = fullfile(args.mask_dir, sprintf('%s.nii',maskName));

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

cd(args.bold_dir)

Preview_raw_filenames = {'Preview1_corr_dt_mcf_brain.nii','Preview2_corr_dt_mcf_brain.nii'};
DFencode_raw_filenames={'DF_encoding_1_corr_dt_mcf_brain.nii', 'DF_encoding_2_corr_dt_mcf_brain.nii'};

subj = load_spm_pattern(subj, 'Preview_epis', maskName, Preview_raw_filenames);
subj = load_spm_pattern(subj, 'DFencode_epis', maskName, DFencode_raw_filenames);

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

for x = 1:length(rsa.DFencode.trialnum);
    if x >   args.DFencode.nTRs/2
        all_DFencode_runs_selector(x) = 2;
    else all_DFencode_runs_selector(x) = 1;
    end
end

subj = initset_object(subj, 'selector', 'all_DFencode_runs', all_DFencode_runs_selector);

  fprintf('[+] z-scoring DFencode data\n\n');
  subj = zscore_runs(subj,'DFencode_epis', 'all_DFencode_runs');
  

%*************** option: read epis in mask | wholebrain - NOT DONE
%currently everything is done through mask
 %*************** define: args.wholebrain '0 | 1' - NOT DONE
 
 
%% ============= Shift regressors and average pattern


%% ============== Average Patterns for option 1 (rsa_type)
%  if rsa_type = 1
%      
%      % take regressors and use training data to shift 2-3 TRs over
%  else
     

%identifying the first TR of each trial for preview
for k = 1:args.trialnum
    trialTR_idx = find(rsa.preview.trialnum == k);
    %for shift
    trialTR_idx_sh = trialTR_idx+args.shiftTR;
    trialTR_mean = trialTR_idx_sh(args.preview.meanTR_start:(args.preview.meanTR_start+args.preview.meanTR_length-1));
    trial_patt = mean(subj.patterns{1,3}.mat(:,trialTR_mean(1):(trialTR_mean(args.preview.meanTR_length))),2);
    rsa.preview.mean.patterns(:,k) = trial_patt;
end



    for m = 1:args.trialnum
%         trialTR_idx = find(rsa.DFencode.trialnum == m);
%         %for shift
%         trialTR_idx_sh = trialTR_idx+args.shiftTR;
%         trialTR_mean = trialTR_idx_sh(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
%         trial_patt = mean(subj.patterns{1,4}.mat(:,trialTR_mean(1):(trialTR_mean(args.DFencode.meanTR_length))),2);
%         rsa.DFencode.mean.patterns(:,m) = trial_patt;
        
            trialTR_idx_match = find(rsa.DFencode.DFencode2preview == m);
            trialTR_idx_match_sh = trialTR_idx_match+args.shiftTR;
            trialTR_mean_match = trialTR_idx_match_sh(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
            trial_patt_match = mean(subj.patterns{1,4}.mat(:,trialTR_mean_match(1):(trialTR_mean_match(args.DFencode.meanTR_length))),2);
           rsa.DFencode.mean.patterns_match(:,rsa.DFencode.DFencode2preview_nonexpanded(m)) = trial_patt_match;
        
        %     trialTR_idx_match = find(rsa.DFencode.DFencode2preview == m);
        %     trialTR_mean_match = trialTR_idx_match(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
        %     trial_patt_match = mean(subj.patterns{1,4}.mat(:,trialTR_mean_match(1):(trialTR_mean_match(args.DFencode.meanTR_length))),2);
        %    rsa.DFencode.mean.patterns_match(:,rsa.DFencode.DFencode2preview_nonexpanded(m)) = trial_patt_match;
    end
    
    
    % corr_matrix_match = zeros(args.trialnum, args.trialnum);
    corr_matrix_match_r1 = zeros(args.trialnum/2,args.trialnum/2);
    for x = 1:args.trialnum/2 %trial number
        for y = 1:args.trialnum/2 %trial number
            corr_matrix_match_r1(x, y) = corr2(rsa.preview.mean.patterns(:,x), rsa.DFencode.mean.patterns_match(:,y));
            rsa.DFencode.mean.smatrix.corr_matrix_match_r1 = corr_matrix_match_r1;
        end
    end
    
    %run 2
    corr_matrix_match_r2 = zeros(args.trialnum/2,args.trialnum/2);
    for x = 1:args.trialnum/2 %trial number
        a = x+args.trialnum/2;
        for y = 1:args.trialnum/2 %trial number
            b = y + args.trialnum/2;
            corr_matrix_match_r2(x, y) = corr2(rsa.preview.mean.patterns(:,a), rsa.DFencode.mean.patterns_match(:,b));
            rsa.DFencode.mean.smatrix.corr_matrix_match_r2 = corr_matrix_match_r2;
        end
    end
    
    %full
    
    corr_matrix_match_full = zeros(args.trialnum,args.trialnum);
    for x = 1:args.trialnum %trial number
        
        for y = 1:args.trialnum %trial number
            
            corr_matrix_match_full(x, y) = corr2(rsa.preview.mean.patterns(:,x), rsa.DFencode.mean.patterns_match(:,y));
            rsa.DFencode.mean.smatrix.corr_matrix_match_full = corr_matrix_match_full;
        end
    end
    
    cd(args.output_dir)
    run1_fig = figure;
    set(run1_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_r1); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    
    saveas(run1_fig, sprintf('emodif%s_run1',subjNum),'png')
    
    run2_fig = figure;
    set(run2_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_r2); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(run2_fig, sprintf('emodif%s_run2',subjNum),'png')
    
    full_fig = figure;
    set(full_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_full); colormap('jet'); colorbar;
    
    
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(full_fig, sprintf('emodif%s_runfull',subjNum),'png')



    %for each DFencode TR
    for m = 1:args.trialnum
        
        trialTR_idx_match = find(rsa.DFencode.DFencode2preview == m);
%         for shift
        trialTR_idx_match_sh = trialTR_idx_match+args.shiftTR;
        for trial_TR = 1:length(trialTR_idx_match_sh)
            TR_patt_match = subj.patterns{1,4}.mat(:,trialTR_idx_match_sh(trial_TR));
            

%         for trial_TR = 1:length(trialTR_idx_match)
%             TR_patt_match = subj.patterns{1,4}.mat(:,trialTR_idx_match(trial_TR));
            
            
            rsa.DFencode.trial(m).patterns_match(:,trial_TR) = TR_patt_match;
        end
    end
    
    %for correlation
        
    
    for trial_TR = 1:args.DFencode.trial_length
        for x = 1:args.trialnum %trial number
            
            for y = 1:args.trialnum %trial number
                
                rsa.DFencode.bytrialTR.smatrixbytr(trial_TR).corr_matrix_match(x,y) = corr2(rsa.preview.mean.patterns(:,x), rsa.DFencode.trial(y).patterns_match(:,trial_TR));
            end
        end
    end
    

    % elseif rsa_type == 2
    %
    % end
    
    % heatmap2 = zeros(args.trialnum, args.trialnum);
    % for x = 1:args.trialnum %trial number
    %     for y = 1:args.trialnum %trial number
    %         heatmap2(x, y) = corr2(rsa.preview.mean.patterns(:,x), rsa.DFencode.mean.patterns(:,y));
    %     end
    % end
    %
    % heatmap_fig2 = figure;
    % set(heatmap_fig2, 'Position', [0 0 1500 500])
    %
    % subplot(1,2,1)
    % imagesc(heatmap2); colormap('jet'); colorbar;
    % hold on
    
    
    
    
    rsa.parameters = args;
    results_file = sprintf('%s/%s_%dto%d_rsa_results.mat',...
           args.output_dir, args.subjID, args.DFencode.meanTR_start, args.DFencode.meanTR_length);
       save(results_file,'rsa')
       
    cd (start_dir); 
end
    
    
    
    
    
