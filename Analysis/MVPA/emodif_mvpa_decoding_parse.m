function emodif_mvpa_decoding_parse(subjNum,maskName,test_phase, test_date);
%*emodif_mvpa_decoding_study('101','tempoccfusi_pHg_LOC_combined_epi_space','L2logreg','fsoner', '50','02', 'DFencode', 2, 'yes')
  version = '2018Apr30';
  args.nTRs = 426; % 135 trials per run, 5 TRs per rest, each post miniblock. each single run 145 stim TRs and 75 rest TRs.  
  args.study_nTRs = 426;
  args.preview_nTRs = 366;
 args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  
    %for astoria
%    args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
   args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);

  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s',args.subj_dir, test_phase, test_date);

cd(args.output_dir)
if strcomp(test_phase, 'DF') == true
load(sprintf('emodif_%s_DFencode_class_perf.txt',subjNum));
load(sprintf('emodif_%s_DFencode_regressors.txt',subjNum));
elseif strcomp(test_phase, 'preview') == true
    load(sprintf('emodif_%s_preview_class_perf.txt',subjNum));
    load(sprintf('emodif_%s_preview_regressors.txt',subjNum));
   
end
end