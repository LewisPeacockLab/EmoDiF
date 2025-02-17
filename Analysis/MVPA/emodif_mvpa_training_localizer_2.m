% function [localizer] = imdif_mvpa_training_localizer_2(subjNum,maskName,featSel,fsThresh,classifier,categories,penalty,shiftTRs)
function emodif_mvpa_training_localizer_2(subjNum,maskName,classifier,categories,penalty,shiftTRs,rest_shift)
  
  %----------------------------------------------------------------------
  % [stuff] = emodif_mvpa_training_localizer(... ALL ARGS ARE STRINGS ...)
  % * for development -
  % emodif_mvpa_training_localizer_2('104','tempoccfusi_pHg_LOC_combined_epi_space','L2logreg','fsowr', '50','02', 2)
  % * subjNum     = subject ID (e.g., '110915')
  % * maskName    = name of mask to use to read in data (no SUBJID)
  % * featSel     = 1|0: do voxelwise ANOVA feature selection, p=0.05
  % * fsThresh    = p-value threshold for feature selectin (default 0.05)
  % * classifier  = name of classifier to use ('L2logreg','logreg', or 'ridge')
  % * categories  = 'fsor' --OR choose one short_cond ('f', 's', 'o', or 'r'! (face|scene|object|rest
  % * penalty    = (e.g., 50|75|100) for 'logreg' & 'ridge' only
  % * shiftTRs    = # of TRs to shift data to adjust for hemodynamic lag
  %
  % 1. read in epis.
  % 2. train a classifier on and do cross-validation.
  %
  %----------------------------------------------------------------------
  version = '2018Apr24';
  
  %% setup parameters
  %-----------------------------------------------------------------------%
  % 
  
  % argument checking
  verify_string(subjNum);
  verify_string(maskName);
  %   verify_string(featSel);
  %   verify_string(fsThresh);
  verify_string(classifier);
  verify_string(categories);
%   verify_string(train_phase);
  verify_string(penalty);
  verify_string(shiftTRs);
  rest_shift_string = num2str(rest_shift);
  
  % setup args structure to keep arguments
  clear args;
  [s,user] = system('whoami');
  [s,host] = system('hostname');
   
  args.script = mfilename;
  args.phase = 'localizer';
  args.nTRs = 342; % 108 trials per run, 5 TRs per rest, each post miniblock. each single run 108 stim TRs and 60 rest TRs.  
  
  args.user = strtrim(user); 
  args.hostname = strtrim(host);
  args.datestamp = sprintf('%s_%s',date,args.hostname);
  args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  
  args.maskName = sprintf('%s.nii',maskName);
  %args.maskName = sprintf('%s/%s/BOLD/%s','/corral-repl/utexas/lewpealab', args.subjID, 'mvpa/mask');
  %   args.featSel = str2num(featSel);
  %   args.featSelThresh = str2num(fsThresh);
  
  args.classifier = classifier;
  args.categories = categories;
%   args.trainphase = train_phase;
  args.penalty = str2num(penalty);
  args.shiftTRs = str2num(shiftTRs);
  args.condNames = {'face','scene','object','word' 'rest'};
  args.condNames_short = {'f','s','o','w','r'};
  args.impmapType = 'mcduff';
  
%   for astoria
   args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
%    args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
%for box synch
% basedir = '/Users/tw24955/Box Sync/';
% args.subj_dir = sprintf('%sEmoDF_LewPeaLab/emodif_data/%s', basedir,args.subjID);
   
  
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s/%s',args.subj_dir, args.phase, maskName, date);
  mkdir(args.output_dir);   


  %----------------------------------------------------------------------
  % turn on diary to capture analysis output
  %
%   diary on;
%   diary(sprintf('%s/%s_diary.txt',args.output_dir,args.subjNum));
%   fprintf('###########################################\n\n');
%   disp(args);
%   fprintf('###########################################\n');


%% initialize subject structure with 'study' and 'subject' info
  %-----------------------------------------------------------------------%
  %
  subj = init_subj('emodif', args.subjID);
  
  %% read in the regressors & selectors
  %-----------------------------------------------------------------------%
  % 
  
  % read in subject-specific regressors
  start_dir = pwd;
  cd(args.regs_dir);
  the_regressors = 'EmoDiF_mvpa_allregs.mat';
  %the_regressors = sprintf('/work/03034/twang04/Dropbox/studydata/imdif/data/sub_%s/IMDIF_MVPA_REGS.mat', args.subjNum);
  
  load(the_regressors);
  
  loc_conds = mvpa_regs.localizer.cat;
  loc_runs = mvpa_regs.localizer.run;
  loc_trial = mvpa_regs.localizer.trial;
%   loc_probe = mvpa_regs.localizer.probe;
  %loc_conds : 1|2|3|4 (object,face,scene,rest)
  %loc_runs  : 1|2|3|4 (run # of the trial)
  %loc_trial : 1|2|...|24 (trial # within the run)
  %loc_probe : 1|2 (absent|present)
  %loc_acc   : 0|1 (WM accuracy of the trial)
  %loc_rt    : 0|1 (WM RT of the trial)
        
  %% MANUALLY adjust the regressors & selectors
  %-----------------------------------------------------------------------%
  % 
  
  % remove any runs?
  switch args.subjNum
    case '20140429'
      loc_runs(loc_runs==3) = NaN;
%     case '05_14067'
%       loc_runs(loc_runs==4) = NaN;
  end
  
  % remove any empty TRs
  norun_trs = isnan(loc_runs);
  regs = {'loc_conds', 'loc_runs', 'loc_trial'};
  if ~isempty(norun_trs)
    fprintf('\n\n## removing %d TRs from localizer regressors\n', count(norun_trs));
    for i = 1:length(regs)
      cmd = sprintf('%s(norun_trs) = [];',regs{i});
      eval(cmd);
    end
  end
  
  %% Parse the arguments to determine which categories to use for MVPA
  %-----------------------------------------------------------------------%
  %
  conds_to_use = zeros(1,length(args.condNames));
  for i = 1:size(args.condNames,2)
    if any(ismember(args.categories,args.condNames_short{i}))
      cond = 1;
    else
      cond = 0;
    end
    conds_to_use(i) = cond;
  end
    
  %% Judy built out regressors already %%
  
  my_conds = mvpa_regs.localizer.cat; 
  all_runs = mvpa_regs.localizer.run;
  all_trials = mvpa_regs.localizer.trial;
  all_TRs = mvpa_regs.localizer.TR;
  
  

  
  

  
  
% %   %% build the regressors and selectors and add to subj structure and
% %   evenly sample rest
%   %build an ID for all rest blocks
% 
%   rest_ID_idx = find(my_conds == 6); %should be 156
%   
%   %needs to treat each run separately so pre 213 and post 214 would be the
%   %end of the first block and start of the next block. 
%   
%  rest_ID_idx_run1 = rest_ID_idx(1,1:78);
%  rest_ID_idx_run2 = rest_ID_idx(1,79:end); 
%  
%  %generate 78 random numbers
%  rand1 = rand(1,78);
%  rand2 = rand(1,78);
%  %sort these 78 random numbers with your idx run within run
%  rest_ID_idx_run1_rand = vertcat(rand1, rest_ID_idx_run1);
%  rest_ID_idx_run2_rand = vertcat(rand2, rest_ID_idx_run2);
%  
%  rest_ID_idx_run1_rand = rest_ID_idx_run1_rand';
%  rest_ID_idx_run2_rand = rest_ID_idx_run2_rand';
%  
% rest_ID_idx_run1_rand_sorted = sortrows(rest_ID_idx_run1_rand);
% rest_ID_idx_run2_rand_sorted = sortrows(rest_ID_idx_run2_rand);
% 
%  face_ID_idx = find(my_conds == 1);
%  num_face_ID_idx_perrun = (length(face_ID_idx))/2;
%  
%  num_rest_TR_elim = length(rest_ID_idx_run1) - num_face_ID_idx_perrun;
%  
%       for x = 1:num_rest_TR_elim
%          my_conds(1,rest_ID_idx_run1(x)) = 0;
%          my_conds(1,rest_ID_idx_run2(x)) = 0;
%      end
 

 
%   %-----------------------------------------------------------------------%
%   % 
%  
  all_conds = []; % 5 conditions face, scene, object, word, rest  
%   all_runs = [];  % 1266 - 3 runs
%   all_trials = []; % 1266 - 432 trials
%   all_TRs = []; %1266 - 1-26 of miniblocks.
%   
  %%%

  
  
  
  %build out conditions
  for k = 1:length(conds_to_use) % 
      temp_conds= [];
      for i = 1:length(my_conds)
          if conds_to_use(k) == 1;
              if my_conds(i) == k
                  temp_conds(i)=1;
              else
                  temp_conds(i)=0;
              end
          end
          
      end
      all_conds = vertcat(all_conds,temp_conds);
  end

  %sanity check
for k = 1:sum(conds_to_use)
    count_conds(k)=sum(all_conds(k,:));
end
    

%  % IF rest is used
 rest_num = strfind(categories,'r');
 if conds_to_use(5) == 1 %6 is REST
     
     
     
       %NONRANDOMLY select rest - take  consistent 2 TRs of every rest block and
  %one TR for the second to last rest block. This can shift Rest sampling
  %until we find the best Rest Sample. 27 TRs. 
  
  n_trial_length = 9;
  n_rest_length = 5;
  n_break_length = 8;
  n_rest_block = 11;
  last_rest_block = 12;
  rest_sample = 2;
  block_size = 171;
  
  % 12 rest blocks of 5TRs following 9 trial TRs, 9 trial TR has 8 rest
  % blocks that follow. 
  
  %rest will always sample two TRs, so a shift of 0 will give TRs 1 and 2
  %of rest.  Shift of 1 will give TRs 2 and 3 or rest, Shift of 2 will give
  %TRs 3 and 4 of rest, shift of 3 will give TRs 4 and 5 of rest. Shift
  %maximum is 3.  will sample last 5 TRs following last trial block. 
  
    rest_vector = zeros(1,block_size);
    rest_vector_dummy = zeros(1,block_size);
    rest_vector_new = zeros(1,block_size);

      for i = 1:n_rest_block
          rest_vector(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(1+rest_shift))))=1;
          rest_vector(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(2+rest_shift))))=1;
          
      end
      
      for i = 1:n_rest_block
          rest_vector_new(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(1+rest_shift))))=1;
          rest_vector_new(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(2+rest_shift))))=1;
          
      end
      
      for i = 1:n_rest_block
          rest_vector_dummy(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(1+rest_shift)))+1)=1;
          rest_vector_dummy(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(2+rest_shift)))+1)=1;
          
      end
      
      if  subjNum == '113'
          rest_vector_new(1,5:end) = rest_vector_dummy(1,5:end);
          
      elseif subjNum == '115'
          rest_vector_new(1,102:end) = rest_vector_dummy(1,102:end); %102 so minus 1
          
      elseif subjNum == '117'
          rest_vector_new(1,119:end) = rest_vector_dummy(1,119:end);
          
      elseif subjNum == '118'
          rest_vector_new(1,132:end) = rest_vector_dummy(1,132:end);
          
      elseif subjNum == '124'
          rest_vector_new(1,57:end) = rest_vector_dummy(1,57:end);
   
          
      end
      %         rest_vector(1,((n_trial_length*last_rest_block)+(((n_rest_length*last_rest_block)-n_rest_length)+(1+rest_shift))))=1;
      rest_vector(1,(length(rest_vector)-7):(length(rest_vector)-3))= 1; % setting the last rest vectors. 
      rest_vector_new(1,((length(rest_vector_new)-6):(length(rest_vector_new)-2))) = 1;
      
      if  subjNum == '113'
          new_rest = horzcat(rest_vector_new, rest_vector);
          
      elseif subjNum == '115'
          new_rest = horzcat(rest_vector, rest_vector_new);
          
      elseif subjNum == '117'
          new_rest = horzcat(rest_vector_new, rest_vector);
          
      elseif subjNum == '118'
          new_rest = horzcat(rest_vector, rest_vector_new);
          
      elseif subjNum == '124'
          new_rest = horzcat(rest_vector_new, rest_vector);
      else
          
        new_rest = horzcat(rest_vector, rest_vector);
          
      end

      
        all_conds(rest_num,:)=new_rest;

  end
 
  
  fprintf('\n\n## using stimulus categories: %s\n',args.categories);
  whos all_regressors
  
  % add regressors & selectors to structure
  conds_name = sprintf('%s_conds',args.phase);
  runs_name = sprintf('%s_runs',args.phase);
  trials_name = sprintf('%s_trials',args.phase);
  TRs_name = sprintf('%s_TRs',args.phase);
  
  subj = initset_object(subj,'regressors',conds_name,all_conds);
  subj = initset_object(subj,'selector',runs_name,all_runs);
  subj = initset_object(subj,'regressors',trials_name, all_trials);
  subj = initset_object(subj,'regressors',TRs_name, all_TRs);
  
  %% shift regressors to account for hemodynamic lag
  %-----------------------------------------------------------------------%
  %
  
  subj = shift_regressors(subj, conds_name, runs_name, args.shiftTRs);
  sh_conds_name = sprintf('%s_sh%d',conds_name,args.shiftTRs);
  subj = shift_regressors(subj, trials_name, runs_name, args.shiftTRs);
  sh_trials_name = sprintf('%s_sh%d',trials_name,args.shiftTRs);
  subj = shift_regressors(subj, TRs_name,runs_name, args.shiftTRs);
  sh_TRs_name = sprintf('%s_sh%d',TRs_name,args.shiftTRs);
  
  
  %shift regressors
  shifted_regs = get_mat(subj,'regressors',sh_conds_name);
  shifted_trials = get_mat(subj, 'regressors', sh_trials_name);
  shifted_TRs = get_mat(subj, 'regressors', sh_TRs_name);
  subj = set_mat(subj, 'regressors', sh_conds_name, shifted_regs);
  fprintf('\n\n## # of regressors in each condition: %s\n\n',mat2str(sum(shifted_regs,2)'));
 
%% exclude rest timepoints
  %-----------------------------------------------------------------------%
  %
  % need to exclude rest timepoints from the analysis
  % want to set all 'baseline' timepoints to '0' so they are excluded
  % create a new selector called 'no_rest' with rests points set to 0
  % to exclude not used rest timepoints during testing*
  
  subj = create_norest_sel(subj,sh_conds_name);
  
  norest_sels_name = [sh_conds_name '_norest'];

  
  
  %% read in voxels from the pre-processed files.
  %-----------------------------------------------------------------------%
  %
  
  cd(args.mask_dir);
  
  fprintf('\n\n## reading in %s data from %s\n',args.phase, args.maskName);
  subj = load_spm_mask(subj,'read_mask',args.maskName);
  
  cd (args.bold_dir);
  
  % determine how many epis there are
  % only grab data from runs that are selected
  runs_to_use = unique(loc_runs);
%   cmd = 'ls -1 imdif_localizer*/bold_dt_mcf_brain.nii | grep -c local'; %dt is 30 sigma dt128 is 128 sigma 
%   [s,r] = system(cmd); %using the system to get this
%   nEPIs = str2num(r);
%   assert(nEPIs == count(runs_to_use),'** Mismatch on # of runs to use');
% 
%   
%   raw_filenames={};
%   cmd = 'ls -1 imdif_localizer_epi_*/bold_dt_mcf_brain.nii | cut -d_ -f4 | cut -d/ -f1 | sort -nr';
%   [~,r] = system(cmd);
%   %listing files
%   r = strtrim(r);
%   raw_filenum= strsplit(r);
%   for i = 1:length(raw_filenum);
%   raw_filename = sprintf('imdif_localizer_epi_%s/bold_dt_mcf_brain.nii', raw_filenum{i});
%   raw_filenames = [raw_filename raw_filenames];
%   end
%       
%       
%   raw_filenames=raw_filenames(1:count(runs_to_use));

  raw_filenames = {'MVPA_training_1_corr_dt_mcf_brain.nii', 'MVPA_training_2_corr_dt_mcf_brain.nii'};
  
  readmask = get_mat(subj,'mask','read_mask');
  
  % read in EPI data as 'doubles'
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.phase, args.maskName, count(readmask));
  subj = load_spm_pattern(subj, 'localizer_epis', 'read_mask', raw_filenames);
  
  
  % Special subject handling: DATA
  
  %PAD Pattern here if missing volumes (e.g. 04_140626 is missing 60
  %volumes from the first localizer task) 1536 is hardcoded as the number
  %of volumes
  
  % %   switch args.subjNum
  % %     case '04_140626'
  % %       % Missing fMRI data in run 1 TRs 325-384.
  % %       % Pad data with NaNs for missing TRs.
  % %
  % %       total_TRs = 1536;
  % %
  % %       epis = get_mat(subj, 'pattern', 'localizer_epis');
  % %       [n_vox,n_trs] = size(epis);
  % %       assert(n_trs==1476); % 60 TRs short in run 1
  % %
  % %       nanpad = NaN(n_vox,(total_TRs-n_trs));
  % %       new_epis = [epis(:,1:324) nanpad(:,1:60) epis(:,325:end)];
  % %       subj = set_mat(subj, 'localizer_epis', new_epis);
  % %
  % %   end
  %% z-score data
  %-----------------------------------------------------------------------%
  %
  
  fprintf('\n\n## z-scoring data\n');
  subj = zscore_runs(subj, 'localizer_epis', runs_name);
  
  
  %% create cross-validation indices
  %-----------------------------------------------------------------------%
  %
  % We are going to create a group of selectors in anticipation of the
  % n-minus-one cross-validation scheme that we will use to train and test our
  % classifier.
  
  subj = create_xvalid_indices(subj,runs_name,'actives_selname',norest_sels_name);

  %% feature select anova  TW:not using ANOVA to feature select as of yet.  Using BET mask as feature selection mask. 
  %-----------------------------------------------------------------------%
  %
  xval_mask = 'read_mask';
%   if args.featSel
%     
%     fprintf('\n\n## starting feature selection ANOVA on data\n')
%     subj = feature_select(subj,'localizer_epis_z',sh_conds_name, ...
%       [runs_name '_xval'], 'thresh', args.featSelThresh);
%     xval_mask = ['localizer_epis_z' '_thresh0' fsThresh];
%     
%     % Write feature-selection mask to disk.
%     sample_filename = mask_fullname;
%     
%     for i = 1:length(runs_to_use)
%       mask_name = sprintf('%s_%d',xval_mask,i);
%       output_name = sprintf('%s_%s_tr%d_%s_%d_%d+orig.BRIK',...
%         args.maskName,xval_mask,args.shiftTRs,args.categories,length(runs_to_use),runs_to_use(i));
%       write_to_afni(subj,'mask',mask_name,sample_filename,...
%         'output_filename',output_name,...
%         'overwrite_if_exist',true);
%       delete allvols*
%     end
%     
%   end
  
  %% do cross-validation on the data!
  %-----------------------------------------------------------------------%
  %
  % determine which classifier algorithm we're going to use
  
  switch(args.classifier)
        
    case 'L2logreg'
      
      class_args.train_funct_name = 'train_L2_RLR';
      class_args.test_funct_name = 'test_L2_RLR';
            class_args.penalty = args.penalty;
            class_args.lambda = 'crossvalidation';

%       class_args.lambda = args.penalty;      
    case 'logreg'
      class_args.train_funct_name = 'train_logreg';
      class_args.test_funct_name = 'test_logreg';
      class_args.penalty = penalty;
      
    case 'ridge'
      class_args.train_funct_name = 'train_ridge';
      class_args.test_funct_name = 'test_ridge';
      class_args.penalty = str2double(penalty);
          
    otherwise
      disp('(-) unknown classifier type!');
      return
  
  end
  
  fprintf('\n\n## %s: Running CROSS-VALIDATION\n', mfilename);
    
  [subj results] =  cross_validation(subj, ...
      'localizer_epis_z', sh_conds_name, [runs_name '_xval'], ...
      xval_mask, class_args);
  
    results.args = args;
    % save results structure 
    mkdir(args.output_dir)
    results_file = sprintf('%s/%s_%s_results.mat',...
      args.output_dir, args.subjID, args.phase);
    save(results_file, 'results')
    
    %write output to xls
    concat_acts=[];
    concat_guesses=[];
    concat_desireds=[];
    concat_accuracy=[];
    
    for i=1:length(runs_to_use)
        results_acts = results.iterations(i).acts;
        concat_acts = [concat_acts results_acts];
        results_guesses = results.iterations(i).perfmet.guesses;
        concat_guesses = [concat_guesses results_guesses];
        results_desireds = results.iterations(i).perfmet.desireds;
        concat_desireds = [concat_desireds results_desireds] ;
        results_accuracy = results.iterations(i).perfmet.corrects;
        concat_accuracy = [concat_accuracy results_accuracy] ;
    end
   
    concat_perf=[concat_guesses; concat_desireds; concat_accuracy; concat_acts];
    shifted_trial_TR_reg=[shifted_trials; shifted_TRs; shifted_regs];
    
    dlmwrite(sprintf('%s/%s_%s_concat_perf.txt', args.output_dir, args.subjID, args.phase),concat_perf);
    dlmwrite(sprintf('%s/%s_%s_shifted_info.txt', args.output_dir, args.subjID, args.phase),shifted_trial_TR_reg);
    
    fn = sprintf('%s/%s_%s_parameters.txt',  args.output_dir, args.subjID, args.phase);
    fid=fopen(fn,'w');
    fprintf(fid,sprintf('%s %s %s %s %s %s %s %s %s %s',args.subjID, args.phase, maskName,  classifier, categories, penalty, shiftTRs, rest_shift_string));
%     [rows cols] = size (concat_perf);
%     x= repmat('%.4f\t',1,(cols-1));
%     fprintf(fid,[x,'%.4f\n'],concat_perf);
    fclose(fid);
    
%     subjNum,maskName,classifier,categories,train_phase,penalty,shiftTRs)
%     dlmwrite(sprintf('%s/%s_%s_concat_desireds.txt',args.output_dir, args.subjID, args.phase),concat_desireds)
%     dlmwrite(sprintf('%s/%s_%s_concat_accuracy.txt',args.output_dir, args.subjID, args.phase),concat_accuracy)
%     dlmwrite(sprintf('%s/%s_%s_shifted_regs.txt',args.output_dir, args.subjID, args.phase),shifted_regs)
%     dlmwrite(sprintf('%s/%s_%s_shifted_trials.txt',args.output_dir, args.subjID, args.phase),shifted_trials)
%     dlmwrite(sprintf('%s/%s_%s_shifted_TRs.txt',args.output_dir, args.subjID, args.phase),shifted_TRs)
    
%  GET EXCEL SERVER WORKING and we can export as sheets.  until then, CSV files
%     xlsfile = sprintf('%s/%s_%s_results.xls',...
%         args.output_dir, args.subjID, args.phase);
%     xlsrange = 'A1';
%    sheet ={'cond_concat','shifted_regs', 'shifted_trials', 'shifted_TRs'};
%     data = {concat_acts, shifted_regs, shifted_trials, shifted_TRs};
%  
%     for i = 1:length(data)
%         current_sheet = [sheet{i}];
%         current_data = [data{i}];
%         xlswrite(xlsfile, current_data, current_sheet, xlsrange);
%     end
    
clear results

diary off;	
cd (start_dir); 
end
% 
function [] = verify_string(thisarg)  
  assert(ischar(thisarg),'** %s must be a string', num2str(thisarg));
end

function [scratchpad] = train_L2_RLR(trainpats,traintargs,class_args,cv_args)
 
%% Train a K class logistic regression classifier, with optional regularization via L2 norm of weight vector(s)
 
% To implement L2 regularization, user must specify class_args.lambda parameter
 
 
 
%% process arguments
                          
% default parameters
defaults.regularization = 'L2';
defaults.lambda         = 1;
defaults.optimization   = 'minimize';
defaults.labelsGroup    = []; 
 
 
 class_args = mergestructs(class_args, defaults);
 
  regularization = class_args.regularization;
  lambda         = class_args.lambda;
  optimization   = class_args.optimization;
  lambda         = class_args.penalty;
   
%% call the classifier function
 
model = classifierLogisticRegression( trainpats', traintargs', 'lambda', lambda);
 
%% pack the results
 
scratchpad.W = model.W;
scratchpad.weights = model.weights;
scratchpad.biases  = model.biases;
 
scratchpad.lambda         = lambda;
scratchpad.regularization = regularization;
scratchpad.optimization   = optimization;
end
