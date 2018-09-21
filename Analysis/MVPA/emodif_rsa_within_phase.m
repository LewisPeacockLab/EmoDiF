function [rsa] = emodif_rsa_within_phase(subjNum,maskName, result_date, shift, dfencodeTRstart, dfencodeTRlength)
%  [rsa] = emodif_rsa_within_phase('103','tempoccfusi_pHg_LOC_combined_epi_space', '29-Aug-2018', 2, 2, 3)
% [rsa] = emodif_rsa_within_phase('101','scene-object_05', 'preview', '29-Aug-2018', 2, 2, 3)
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
  args.Localizer.meanTR_start = 4;
  args.Localizer.meanTR_length = 3;
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
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.DFencode_dir = sprintf('%s/results/%s/%s',args.subj_dir, maskName, result_date);
  args.preview_dir = sprintf('%s/results/%s/%s',args.subj_dir, maskName, result_date);
  args.output_dir = sprintf('%s/results/rsa_results/intraphase/%s/%d/%s',args.subj_dir,maskName,args.shiftTR, date);
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

for i = 1:(args.DFencode.trialnum/2)
    x = repmat(i,args.DFencode.trial_length,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for i = ((args.DFencode.trialnum/2)+1):args.DFencode.trialnum
    x = repmat(i,args.DFencode.trial_length,1)';
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


%load masks

mvpa_mask = fullfile(args.mask_dir, sprintf('%s.nii',maskName));

%% ============= Initializing subj. structure:start by creating an empty subj structure
% summarize(subj): summarize all info in subj structure
% get_object/set_object/set_objfield/set_objsubfield
% get_mat/set_mat

subj = init_subj(args.experiment, args.subjID);%identifier of the subj


%% ============= 01: EPI PATTERNS
%*************** load mask + read in epis

subj = load_spm_mask(subj,maskName,mvpa_mask);

mask_voxel = get_mat(subj,'mask', maskName);

fprintf('\n(+) load epi data under mask with %s voxels\n', num2str(count(mask_voxel)));

cd(args.bold_dir)

Preview_raw_filenames = {'Preview1_corr_dt_mcf_brain.nii','Preview2_corr_dt_mcf_brain.nii'};
DFencode_raw_filenames={'DF_encoding_1_corr_dt_mcf_brain.nii', 'DF_encoding_2_corr_dt_mcf_brain.nii'};
Localizer_raw_filenames = {'MVPA_training_1_corr_dt_mcf_brain.nii', 'MVPA_training_2_corr_dt_mcf_brain.nii'};

subj = load_spm_pattern(subj, 'Preview_epis', maskName, Preview_raw_filenames);
subj = load_spm_pattern(subj, 'DFencode_epis', maskName, DFencode_raw_filenames);
subj = load_spm_pattern(subj, 'Localizer_epis', maskName, Localizer_raw_filenames);

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

%% %finding intraphrase patterns. 

%phase patterns
for m = 1:args.preview.trialnum
    preview_trialTR_idx = find(rsa.preview.trialnum == m);
    preview_trialTR_mean = preview_trialTR_idx(args.preview.meanTR_start:(args.preview.meanTR_start+args.preview.meanTR_length-1));
    preview_trial_patt = mean(subj.patterns{1,4}.mat(:,preview_trialTR_mean(1):(preview_trialTR_mean(args.preview.meanTR_length))),2);
    rsa.preview.mean.patterns(:,m) = preview_trial_patt;
end




for m = 1:args.DFencode.trialnum
    DFencode_trialTR_idx = find(rsa.DFencode.trialnum == m);
    DFencode_trialTR_mean = DFencode_trialTR_idx(args.DFencode.meanTR_start:(args.DFencode.meanTR_start+args.DFencode.meanTR_length-1));
    DFencode_trial_patt = mean(subj.patterns{1,5}.mat(:,DFencode_trialTR_mean(1):(DFencode_trialTR_mean(args.DFencode.meanTR_length))),2);
    rsa.DFencode.mean.patterns(:,m) = DFencode_trial_patt;
end



for m = 1:args.Localizer.trialnum
    Localizer_trialTR_idx = find(rsa.Localizer.trialnum == m);
    Localizer_trialCAT = (mvpa_regs.localizer.cat(Localizer_trialTR_idx(1)));
    Localizer_trialTR_mean = Localizer_trialTR_idx(args.Localizer.meanTR_start:(args.Localizer.meanTR_start+args.Localizer.meanTR_length-1));
    Localizer_trial_patt = mean(subj.patterns{1,6}.mat(:,Localizer_trialTR_mean(1):(Localizer_trialTR_mean(args.Localizer.meanTR_length))),2);
    rsa.Localizer.mean.patterns(:,m) = Localizer_trial_patt;
    rsa.Localizer.mean.category(:,m) = Localizer_trialCAT;
end

    
    %full
    
    corr_matrix_match_preview = zeros(args.preview.trialnum,args.preview.trialnum);
    for x = 1:args.preview.trialnum %trial number
        
        for y = 1:args.preview.trialnum %trial number
            
            corr_matrix_match_preview(x, y) = corr2(rsa.preview.mean.patterns(:,x), rsa.preview.mean.patterns(:,y));
            corr_matrix_match_previewz(x,y) = 0.5*log((1+corr_matrix_match_preview(x, y))/(1-corr_matrix_match_preview(x, y)));
            rsa.results.smatrix.corr_matrix_match_preview = corr_matrix_match_preview;
            rsa.results.smatrix.corr_matrix_match_previewz = corr_matrix_match_previewz;
        end
    end
    
    corr_matrix_match_DFencode = zeros(args.DFencode.trialnum,args.DFencode.trialnum);
    for x = 1:args.DFencode.trialnum %trial number
        
        for y = 1:args.DFencode.trialnum %trial number
            
            corr_matrix_match_DFencode(x, y) = corr2(rsa.DFencode.mean.patterns(:,x), rsa.DFencode.mean.patterns(:,y));
            corr_matrix_match_DFencodez(x,y) = 0.5*log((1+corr_matrix_match_DFencode(x, y))/(1-corr_matrix_match_DFencode(x, y)));
            rsa.results.smatrix.corr_matrix_match_DFencode = corr_matrix_match_DFencode;
            rsa.results.smatrix.corr_matrix_match_DFencodez = corr_matrix_match_DFencodez;
        end
    end
    
    corr_matrix_match_Localizer = zeros(args.Localizer.trialnum,args.Localizer.trialnum);
    for x = 1:args.Localizer.trialnum %trial number
        
        for y = 1:args.Localizer.trialnum %trial number
            
            corr_matrix_match_Localizer(x, y) = corr2(rsa.Localizer.mean.patterns(:,x), rsa.Localizer.mean.patterns(:,y));
            corr_matrix_match_Localizerz(x,y) = 0.5*log((1+corr_matrix_match_Localizer(x, y))/(1-corr_matrix_match_Localizer(x, y)));
            rsa.results.smatrix.corr_matrix_match_Localizer = corr_matrix_match_Localizer;
            rsa.results.smatrix.corr_matrix_match_Localizerz = corr_matrix_match_Localizerz;
        end
    end
    
    %%%%% added Localizer - extract out by category - by category
    
    rsa.Localizer.mean.patterns_tosort = vertcat(rsa.Localizer.mean.category,rsa.Localizer.mean.patterns)';
    rsa.Localizer.mean.patterns_sorted = sortrows(rsa.Localizer.mean.patterns_tosort)';
    rsa.Localizer.mean.patterns_bycategory = rsa.Localizer.mean.patterns_sorted(2:end, :);
    
    corr_matrix_match_Localizer_bycat = zeros(args.Localizer.trialnum,args.Localizer.trialnum);
    
        for x = 1:args.Localizer.trialnum %trial number
        
            for y = 1:args.Localizer.trialnum %trial number
                
                corr_matrix_match_Localizer_bycat(x, y) = corr2(rsa.Localizer.mean.patterns_bycategory (:,x), rsa.Localizer.mean.patterns_bycategory (:,y));
                corr_matrix_match_Localizer_bycatz(x,y) = 0.5*log((1+corr_matrix_match_Localizer_bycat(x, y))/(1-corr_matrix_match_Localizer_bycat(x, y)));
                rsa.results.smatrix.corr_matrix_match_Localizer_bycat = corr_matrix_match_Localizer_bycat;
                rsa.results.smatrix.corr_matrix_match_Localizer_bycatz = corr_matrix_match_Localizer_bycatz;
            end
        end
   
    
%     localizer_table = vertcat(rsa.Localizer.mean.category, rsa.results.smatrix.corr_matrix_match_Localizer)';
%     localizer_table_sorted = sortrows(localizer_table)';
%     localizer_corr_matrix_bycat = localizer_table_sorted(2:end, :);
%     vert_cat = rsa.Localizer.mean.category';
%     localizer_corr_matrix_bycat_vert = horzcat(vert_cat, localizer_corr_matrix_bycat);
    
    
    
    
    cd(args.output_dir)
%     run1_fig = figure;
%     set(run1_fig, 'Position', [0 0 1500 1500])
%     
%     subplot(1,1,1)
%     imagesc(corr_matrix_match_r1); colormap('jet'); colorbar;
%     
%     
%     xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
%     ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
%     
%     saveas(run1_fig, sprintf('emodif%s_run1',subjNum),'png')
%     
%     run2_fig = figure;
%     set(run2_fig, 'Position', [0 0 1500 1500])
%     
%     subplot(1,1,1)
%     imagesc(corr_matrix_match_r2); colormap('jet'); colorbar;
%     
%     
%     xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d run 1',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
%     ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d run 1',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
%     saveas(run2_fig, sprintf('emodif%s_run2',subjNum),'png')

    preview_fig = figure;
    set(preview_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_preview); colormap('jet'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Preview raw patterns averaged over TR %d: TR %d ',args.preview.meanTR_start,args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('Preview raw patterns averaged over TR %d: TR %d ',args.preview.meanTR_start, args.preview.meanTR_end),'FontSize',15,'FontWeight','bold');
    
    saveas(preview_fig, sprintf('emodif%s_run_Preview_TRstart%d_mean%d',subjNum,args.preview.meanTR_start,args.preview.meanTR_length),'png')

    DFencode_fig = figure;
    set(DFencode_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_DFencode); colormap('jet'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d ',args.DFencode.meanTR_start,args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d ',args.DFencode.meanTR_start, args.DFencode.meanTR_end),'FontSize',15,'FontWeight','bold');
    
    saveas(DFencode_fig, sprintf('emodif%s_run_DFencode_TRstart%d_mean%d',subjNum,args.DFencode.meanTR_start,args.DFencode.meanTR_length),'png')
    
    Localizer_fig = figure;
    set(Localizer_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(corr_matrix_match_Localizer); colormap('jet'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Localizer raw patterns averaged over TR %d: TR %d ',args.Localizer.meanTR_start,args.Localizer.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('Localizer raw patterns averaged over TR %d: TR %d ',args.Localizer.meanTR_start, args.Localizer.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(Localizer_fig, sprintf('emodif%s_run_Localizer_TRstart%d_mean%d',subjNum,args.Localizer.meanTR_start,args.Localizer.meanTR_length),'png')
    
    
    Localizerbycat_fig = figure;
    set(Localizerbycat_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(rsa.results.smatrix.corr_matrix_match_Localizer_bycat); colormap('jet'); colorbar; 
    args.confplot.origcoloraxis = caxis;
    
    if strcmp(maskName,'JC_Combine_PHc_epi_space') == 1
        caxis([-.25 .25]);
    else
        caxis([-.35 .35]);
    end
    args.confplot.finalcoloraxis = caxis;
    
    ylabel(sprintf('Localizer raw patterns averaged over TR %d: TR %d bycat',args.Localizer.meanTR_start,args.Localizer.meanTR_end),'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('Localizer raw patterns averaged over TR %d: TR %d bycat',args.Localizer.meanTR_start, args.Localizer.meanTR_end),'FontSize',15,'FontWeight','bold');
    saveas(Localizerbycat_fig, sprintf('emodif%s_run_Localizer_bycat_TRstart%d_mean%d',subjNum,args.Localizer.meanTR_start,args.Localizer.meanTR_length),'png')
% 
%     %for each DFencode TR
%     for m = 1:args.trialnum
%         
%         trialTR_idx_match = find(rsa.DFencode.DFencode2preview == m);
% %         for shift
%         trialTR_idx_match_sh = trialTR_idx_match+args.previewshiftTR;
%         for trial_TR = 1:length(trialTR_idx_match_sh)
%             TR_patt_match = subj.patterns{1,4}.mat(:,trialTR_idx_match_sh(trial_TR));
%             
% 
% %         for trial_TR = 1:length(trialTR_idx_match)
% %             TR_patt_match = subj.patterns{1,4}.mat(:,trialTR_idx_match(trial_TR));
%             
%             
%             rsa.DFencode.trial(m).patterns_match(:,trial_TR) = TR_patt_match;
%         end
%     end
%     
%     %for correlation
%         
%     
%     for trial_TR = 1:args.DFencode.trial_length
%         for x = 1:args.trialnum %trial number
%             
%             for y = 1:args.trialnum %trial number
%                 
%                 rsa.DFencode.bytrialTR.smatrixbytr(trial_TR).corr_matrix_match(x,y) = corr2(rsa.preview.mean.patterns(:,x), rsa.DFencode.trial(y).patterns_match(:,trial_TR));
%             end
%         end
%     end
%     

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
    results_file = sprintf('%s/%s_TR%dto%d_rsa_results.mat',...
           args.output_dir, args.subjID, args.DFencode.meanTR_start, args.DFencode.meanTR_end);
       save(results_file,'rsa')
       
    cd (start_dir); 
end
    
    
    
    
    
