function emodif_mvpa_decoding_parse(subjNum,maskName,test_phase, test_date);


  args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s',args.subj_dir, test_phase, test_date);

cd(args.output_dir)
if strcomp(test_phase, 'DF') == true
load(sprintf('emodif_%s_DFencode_class_perf.txt',subjNum);
load(sprintf('emodif_%s_DFencode_regressors.txt',subjNum);
elseif strcomp(test_phase, 'preview') == true
    load(sprintf('emodif_%s_preview_class_perf.txt',subjNum);
    load(sprintf('emodif_%s_preview_regressors.txt',subjNum);
end