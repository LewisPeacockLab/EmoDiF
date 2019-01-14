function [correl] = emodif_rsa_catcompare(maskName)
% this is to compare two MVPA outputs FSOR and FSOWR TR by TR
% emodif_rsa_catcompare('tempoccfusi_pHg_LOC_combined_epi_space')
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
  args.test_phase = 'DFencode';
  args.DFencode.nTRs = 426;
  args.DFencode.trial_length = 7;
  args.DFencode.trial_break = 3;
  args.trialnum = 60;
  
  %for astoria
  args.main_dir = '/Users/tw24955/emodif_data/';
  args.output_dir = sprintf('%s/aggregate_results',args.main_dir);
  mkdir(args.output_dir);
  %for tigger
  %args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
 
  %
%%%% LOAD RESULTFILE %


[x,subj_folders] = xlsread(sprintf('%s/subj_folders.xlsx',args.main_dir)); % just the same for the first 4 subjects

%%%% 
category_comp_correl = [];
groupdata = [];

for s = 1:length(subj_folders)
    
  args.subjID = subj_folders{s,1};
  args.folder_FSOR = subj_folders{s,2};
  args.folder_FSOWR = subj_folders{s,3};
  
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
    %data directories
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  %   args.DFencode_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.test_phase, test_date);
%   args.preview_dir = sprintf('%s/results/%s/%s',args.subj_dir, args.train_phase, train_date);

  args.data_FSOR = sprintf('%s/results/DFencode/%s/%s',args.subj_dir,maskName, args.folder_FSOR);
  args.data_FSOWR = sprintf('%s/results/DFencode/%s/%s',args.subj_dir,maskName, args.folder_FSOWR);

  subjNum = args.subjID;
  
  %load data
  load(sprintf('%s/%s_DFencode_results.mat',args.data_FSOR,subjNum));
  results_FSOR = results;
  clear results;
  load(sprintf('%s/%s_DFencode_results.mat',args.data_FSOWR,subjNum));
  results_FSOWR = results;
  
  category_comp_correl{s,1} = subjNum;
  category_comp_correl{s,2} = corr2(results_FSOR.iterations.acts(1,:),results_FSOWR.iterations.acts(1,:));
  category_comp_correl{s,3} = corr2(results_FSOR.iterations.acts(2,:),results_FSOWR.iterations.acts(2,:));
  category_comp_correl{s,4} = corr2(results_FSOR.iterations.acts(3,:),results_FSOWR.iterations.acts(3,:));
  category_comp_correl{s,5} = corr2(results_FSOR.iterations.acts(4,:),results_FSOWR.iterations.acts(5,:));
  
  groupdata.face{s,1} = results_FSOR.iterations.acts(1,:);
  groupdata.face{s,2} = results_FSOWR.iterations.acts(1,:);
  
  groupdata.scene{s,1} = results_FSOR.iterations.acts(2,:);
  groupdata.scene{s,2} = results_FSOWR.iterations.acts(2,:);
  
  groupdata.object{s,1} = results_FSOR.iterations.acts(3,:);
  groupdata.object{s,2} = results_FSOWR.iterations.acts(3,:);
  
  groupdata.rest{s,1} = results_FSOR.iterations.acts(4,:);
  groupdata.rest{s,2} = results_FSOWR.iterations.acts(5,:);
end


    % hold on
    cat_comp.corr = category_comp_correl;
    cat_comp.parameters = args;
    cat_comp.groupdata = groupdata;
    cat_comp.subj_folders = subj_folders;
    
    results_file = sprintf('%s/comp_compare_results.mat',...
           args.output_dir);
       save(results_file,'cat_comp')
       
    cd (start_dir); 
end
    
    
    
    
    
