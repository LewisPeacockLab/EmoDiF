function [rsa_stuff] = emodif_rsa_preview_local(args, dir)
%RSA set up script based off of hyojeong's clearmem matlab RSA script. ***
%requires princeton toolbox ***

%this RSA is to compare PREVIEW scene related activity with ENCODING scene
%related activity following the same word presentation. 

%Experiment:EmoDF
%Type : Emotional Directed Forgetting with Scene taggs
%Phase: Preview to Study
% Date      : August 2018
% Version   : 1
% Author    : Tracy Wang
% Contact   : tracy.wang@utexas.edu

%argument checking

  verify_string(subjNum);
  
  args.train_phase = 'preview';
  args.test_phase = 'DFencode';
  args.preview.nTRs = 366;
  args.DFencode.nTRs = 426;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  args.preview.trial_length = 6;
  args.study.trial_length = 7;
  
  %for astoria
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
  %for tigger
  args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
  
  %data directories
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/rsa_results/%s/%s',args.subj_dir, args.test_phase, date);
  mkdir(args.output_dir);
  
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

%%% expanding RSA parameters %%%



