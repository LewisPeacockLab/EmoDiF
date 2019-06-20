function emodif_mvpa_localizer_parse(subjNum,phase, maskName, test_date)

%emodif_mvpa_localizer_parse('101','localizer', 'tempoccfusi_pHg_LOC_combined_epi_space', '25-Sep-2018')

  args.break = 3; %how many TRs are at the end of each run regardless of condition
  args.local.nTRs = 426; % 135 trials per run, 5 TRs per rest, each post miniblock. each single run 145 stim TRs and 75 rest TRs.  
  args.local.trialTR = 9;
  
  args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  
  args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
  
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s/%s',args.subj_dir, phase, maskName, test_date);
  args.script_dir = pwd;
  
  
  cd(args.output_dir)
  class_perf = load(sprintf('emodif_%s_localizer_concat_perf.txt',subjNum));
  %class performance is guesses, desireds, accuracy, acts (face scene
  %object word rest)
  regressors = load(sprintf('emodif_%s_localizer_shifted_info.txt',subjNum));
  
  %shifted trials, shifted TRs, shifted regs (which is
  fparameters = fopen(sprintf('emodif_%s_localizer_parameters.txt',subjNum));
  
  %look at word, scene and face across the entire localizer period
  
  desired = class_perf(2,:);
  class_perf_face = class_perf(4,:);
  class_perf_scene = class_perf(5,:);
  class_perf_object= class_perf(6,:);
  class_perf_word = class_perf(7,:);
  class_perf_rest = class_perf(8,:);
  
  
  
  [corr_ws, corrP_ws] = corr(class_perf_word', class_perf_scene');
  [corr_wf, corrP_wf] = corr(class_perf_word', class_perf_face');
  
  %start to parcellate according to miniblock type
  
  face_idx = find(desired' == 1);
  scene_idx = find(desired' == 2);
  object_idx = find(desired' == 3);
  word_idx = find(desired' == 4);
  rest_idx = find(desired' == 5);
  
  
  
  %use desired
  
  %face
  for face_t = 1:length(face_idx)
      
      face_face(face_t,1) = class_perf_face(face_idx(face_t));
      face_scene(face_t,1) = class_perf_scene(face_idx(face_t));
      face_word(face_t,1) = class_perf_word(face_idx(face_t));
      face_object(face_t,1) = class_perf_object(face_idx(face_t));
      face_rest(face_t,1) = class_perf_rest(face_idx(face_t));
      
  end
  
  
  for scene_t = 1:length(scene_idx)
      
      scene_face(scene_t,1) = class_perf_face(scene_idx(scene_t));
      scene_scene(scene_t,1) = class_perf_scene(scene_idx(scene_t));
      scene_word(scene_t,1) = class_perf_word(scene_idx(scene_t));
      scene_object(scene_t,1) = class_perf_object(scene_idx(scene_t));
      scene_rest(scene_t,1) = class_perf_rest(scene_idx(scene_t));
      
  end
  
  for object_t = 1:length(object_idx)
      
      scene_face(object_t,1) = class_perf_face(object_idx(object_t));
      object_scene(object_t,1) = class_perf_scene(object_idx(object_t));
      object_word(object_t,1) = class_perf_word(object_idx(object_t));
      object_object(object_t,1) = class_perf_object(object_idx(object_t));
      object_rest(object_t,1) = class_perf_rest(object_idx(object_t));
      
  end
  
  for word_t = 1:length(word_idx)
      
      word_face(word_t,1) = class_perf_face(word_idx(word_t));
      word_scene(word_t,1) = class_perf_scene(word_idx(word_t));
      word_word(word_t,1) = class_perf_word(word_idx(word_t));
      word_object(word_t,1) = class_perf_object(word_idx(word_t));
      word_rest(word_t,1) = class_perf_rest(word_idx(word_t));
      
  end
  
  for rest_t = 1:length(rest_idx)
      
      rest_face(rest_t,1) = class_perf_face(rest_idx(rest_t));
      rest_scene(rest_t,1) = class_perf_scene(rest_idx(rest_t));
      rest_word(rest_t,1) = class_perf_word(rest_idx(rest_t));
      rest_object(rest_t,1) = class_perf_object(rest_idx(rest_t));
      rest_rest(rest_t,1) = class_perf_rest(rest_idx(rest_t));
      
  end
  
  %correlations
  
  %title reflects blocks - 
  
  [corr_face_ws, corrP_face_ws] = corr(face_word, face_scene);
  [corr_face_wf, corrP_face_wf] = corr(face_word, face_face);
  
  [corr_word_ws, corrP_word_ws] = corr(word_word, word_scene);
  [corr_word_wf, corrP_word_wf] = corr(word_word, word_face);
  
  [corr_scene_ws, corrP_scene_ws] = corr(scene_word, scene_scene);
  [corr_scene_wf, corrP_scene_wf] = corr(scene_word, scene_face);
  
  localizer.data.class_perf = class_perf;
  localizer.corr.corr_face_ws = corr_face_ws;
  localizer.corr.corrP_face_ws = corrP_face_ws;
  localizer.corr.corr_face_wf = corr_face_wf;
  localizer.corr.corrP_face_wf = corrP_face_wf;

  localizer.corr.corr_word_ws = corr_word_ws;
  localizer.corr.corrP_word_ws = corrP_word_ws;
  localizer.corr.corr_word_wf = corr_word_wf;
  localizer.corr.corrP_word_wf = corrP_word_wf;
  
  
  localizer.corr.corr_scene_ws = corr_scene_ws;
  localizer.corr.corrP_scene_ws = corrP_scene_ws;
  localizer.corr.corr_scene_wf = corr_scene_wf;
  localizer.corr.corrP_scene_wf = corrP_scene_wf;
  localizer.corr.corr_all_ws = corr_ws;
  localizer.corr.corrP_all_ws = corrP_ws;
  localizer.corr.corr_all_wf = corr_wf;
  localizer.corr.corrP_all_wf = corrP_wf;
  
  
  
  filename = 'localizer_parsed_corr';
  save(filename, 'localizer');
end
  
  
  
  
  