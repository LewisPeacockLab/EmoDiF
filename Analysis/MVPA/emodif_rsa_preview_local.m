function [rsa_stuff] = emodif_rsa_preview_local(dir)
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
  
  %%% parameters %%%
  args.experiment = 'emodif';
  
  
  %%% end parameters %%%
  
  
  args.train_phase = 'preview';
  args.test_phase = 'DFencode';
  args.preview.nTRs = 366;
  args.DFencode.nTRs = 426;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  args.preview.trial_length = 6;
  args.study.trial_length = 7;
  args.preview.trial_break = 3;
  
  %for astoria
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
  %for tigger
  %args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
  
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

%%% expanding RSA parameters with mvpa_regs %%%

rsa.preview.preview2study = RSA_params.preview2study;


%preview trial numbers

firstblocktrialnum = [];
secondblocktrialnum = [];
for i = 1:30
    x = repmat(i,6,1)';
    firstblocktrialnum = horzcat(firstblocktrialnum, x);
end

for i = 31:60
    x = repmat(i,6,1)';
    secondblocktrialnum = horzcat(secondblocktrialnum, x);
end

previewtrialbreak = zeros(3,1)';

rsa.preview.trialnum = horzcat(firstblocktrialnum,previewtrialbreak,secondblocktrialnum,previewtrialbreak);
% bring DFencode category and responses into preview regressors
for i = 1:args.preview.nTRs
    if rsa.preview.preview2study(i) == 0
        trial_cat = NaN;
    else trial_dfencode = find(mvpa_regs.DFEncode.trial == rsa.preview.preview2study(i));
     trial_cat = mvpa_regs.DFEncode.cat(1,trial_dfencode(1));
     trial_resp = mvpa_regs.DFEncode.subresp(1,trial_dfencode(1));
    end
    rsa.preview.DFencodecat(i) = trial_cat;
    rsa.preview.DFencodeResp(i) = trial_resp;
end

%complete DFencode regressors for RSA
rsa.DFencode.study2preview = RSA_params.study2preview;

for i = 1:args.DFencode.nTRs
    if rsa.DFencode.study2preview(i) == 0
        trial_cat = NaN;
    else trial_study = find(mvpa_regs.preview.trial == rsa.DFencode.study2preview(i));
        trial_cat = mvpa_regs.DFEncode.cat(1,trial_study(1));
        trial_resp = mvpa_regs.DFEncode.subresp(1,trial_study(1));
        
    end
    rsa.DFencode.cat(i) = trial_cat;
    rsa.DFencode.resp(i) = trial_resp;
end
        

%%%%% Making filtered Masks FOR RSA %%%%%%
%apply masks from MVPA - use EPI masks
%load masks

mask_name 


%%% extract PATTERNS %%%%







