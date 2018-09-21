% function [localizer] = imdif_mvpa_training_localizer_tw(subjNum,maskName,featSel,fsThresh,classifier,categories,penalty,shiftTRs)
function emodif_mvpa_decoding_study_noSM_1(subjNum,maskName,classifier,categories,penalty,shiftTRs,test_phase,rest_shift,word_combine)
  
  %----------------------------------------------------------------------
  % [stuff] = emodif_mvpa_training_localizer(... ALL ARGS ARE STRINGS ...)
  % * for development -
  % emodif_mvpa_decoding_study('101','tempoccfusi_pHg_LOC_combined_epi_space','L2logreg','fsoner', '50','02','DFencode', 2, 'yes')
  % * subjNum     = subject ID (e.g., '110915')
  % * maskName    = name of mask to use to read in data (no SUBJID)
  % * featSel     = 1|0: do voxelwise ANOVA feature selection, p=0.05
  % * fsThresh    = p-value threshold for feature selectin (default 0.05)
  % * classifier  = name of classifier to use ('L2logreg','logreg', or 'ridge')
  % * categories  = 'fsor' --OR choose one short_cond ('f', 's', 'o', or 'r'! (face|scene|object|rest
  % * penalty    = (e.g., 50|75|100) for 'logreg' & 'ridge' only
  % * shiftTRs    = # of TRs to shift data to adjust for hemodynamic lag
  % * test_phase = 'DFencode' or 'Preview'
  % * rest_shift = '0-3' out of 5 TRs, we sample 2 for each rest phase
  % * word_combine = combine negative and neutral word evidence 'yes' or
  % 'no'
  % 1. read in epis.
  % 2. train a classifier on and do cross-validation.
  %
  %----------------------------------------------------------------------
  version = '2018Apr30';
  
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
  
  args.train_phase = 'localizer';
  args.test_phase = test_phase;
  [s,user] = system('whoami');
  [s,host] = system('hostname');
   
  args.script = mfilename;

  args.nTRs = 426; % 135 trials per run, 5 TRs per rest, each post miniblock. each single run 145 stim TRs and 75 rest TRs.  
  args.study_nTRs = 426;
  args.preview_nTRs = 366;
  
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
  args.condNames = {'face','scene','object','neutral', 'negative' 'rest'};
  args.condNames_short = {'f','s','o','n','e','r'};
  args.impmapType = 'mcduff';
  
  %for astoria
   args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
%    args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);
   
  
  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s/%s',args.subj_dir, args.test_phase, maskName, date);
  mkdir(args.output_dir);   


  %----------------------------------------------------------------------
  % turn on diary to capture analysis output
  %
  diary on;
  diary(sprintf('%s/%s_diary.txt',args.output_dir,args.subjNum));
  fprintf('###########################################\n\n');
  disp(args);
  fprintf('###########################################\n');


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
  
    %%%%%%%%%%%%%%%%%%%%% SPECIAL CONDITIONS %%%%%%%%%%%%%%%%%
if subjNum == '101'; %101 had an extra TR
    my_conds(218)=6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  

  
  
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
for k = 1:length(conds_to_use)
    count_conds(k)=sum(all_conds(k,:));
end
  

      
    

%  % IF rest is used
 if conds_to_use(6) == 1; %6 is REST
     
     
     
       %NONRANDOMLY select rest - take  consistent 2 TRs of every rest block and
  %one TR for the second to last rest block. This can shift Rest sampling
  %until we find the best Rest Sample. 
  
  n_trial_length = 9;
  n_rest_length = 5;
  n_break_length = 8;
  n_rest_block = 13;
  last_rest_block = 14;
  rest_sample = 2;
  block_size = 213;
  
  % 14 rest blocks of 5TRs following 9 trial TRs, 9 trial TR has 7 rest
  % blocks that follow. 
  
  %rest will always sample two TRs, so a shift of 0 will give TRs 1 and 2
  %of rest.  Shift of 1 will give TRs 2 and 3 or rest, Shift of 2 will give
  %TRs 3 and 4 of rest, shift of 3 will give TRs 4 and 5 of rest. Shift
  %maximum is 3.  
  
  rest_vector = zeros(1,block_size);

  
  for i = 1:n_rest_block
      rest_vector(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(1+rest_shift))))=1;
      rest_vector(1, ((n_trial_length*i)+(((n_rest_length*i)-n_rest_length)+(2+rest_shift))))=1;

  end
        rest_vector(1,((n_trial_length*last_rest_block)+(((n_rest_length*last_rest_block)-n_rest_length)+(1+rest_shift))))=1;
        
        new_rest= horzcat(rest_vector,rest_vector);
        
        all_conds(6,:)=new_rest;
 end
 
   %%% SPECIAL CASES - EXTRA TRS and SUCH ***
  
  if subjNum == '101'
      %TR 218 has nothing in it. 
     all_conds(:,218) = nan;
     my_conds(:,218) = nan;
      
      %this is a special case, 
  end
  
        
 %%%%% COMBINING NEUTRAL AND NEGATIVE WORDS %%%%%
 
 if strcmp(word_combine,'yes') == true; %4 is Words
     args.categories = 'fsowr';
     
     
     rand1 = rand(1,27);
     rand2 = rand(1,27);
     rand3 = rand(1,27);
     rand4 = rand(1,27);
     
     neu_word_ID_idx = find(my_conds == 4); %should be 54
     neg_word_ID_idx = find(my_conds == 5); %should be 54
     
     %needs to treat each run separately so pre 213 and post 214 would be the
     %end of the first block and start of the next block.
     
     neu_ID_idx_run1 = neu_word_ID_idx(1,1:27);
     neu_ID_idx_run2 = neu_word_ID_idx(1,28:end);
     neg_ID_idx_run1 = neg_word_ID_idx(1,1:27);
     neg_ID_idx_run2 = neg_word_ID_idx(1,28:end);


         
         %for the 1st block, we want 14 neutral and 13 negative words
         
         %for the 2nd block, we want 13 neutral and 14 negative words 
         

         neu_ID_idx_run1_rand = vertcat(rand1, neu_ID_idx_run1)';
         neu_ID_idx_run2_rand = vertcat(rand2, neu_ID_idx_run2)';
         
         neg_ID_idx_run1_rand = vertcat(rand3, neg_ID_idx_run1)';
         neg_ID_idx_run2_rand = vertcat(rand4, neg_ID_idx_run2)';
         
         %
         
         %
         neu_ID_idx_run1_rand_sorted = sortrows(neu_ID_idx_run1_rand);
         neu_ID_idx_run2_rand_sorted = sortrows(neu_ID_idx_run2_rand);
         
         neg_ID_idx_run1_rand_sorted = sortrows(neg_ID_idx_run1_rand);
         neg_ID_idx_run2_rand_sorted = sortrows(neg_ID_idx_run2_rand);
         
         % to eliminate 13 neutral (run 1) and 13 negative (run 2) 
         
         for x = 1:13
             all_conds(4,neu_ID_idx_run1_rand_sorted(x,2)) = 0;
             all_conds(5,neg_ID_idx_run2_rand_sorted(x,2)) = 0;
         end
         
         for x = 1:14
             all_conds(4,neu_ID_idx_run2_rand_sorted(x,2)) = 0;
             all_conds(5,neg_ID_idx_run1_rand_sorted(x,2)) = 0;
         end
         
         all_conds_words = (all_conds(4,:)+all_conds(5,:));
        

 %build regressors again
 
 new_all_conds(1,:) = all_conds(1,:);

  new_all_conds(2,:) = all_conds(2,:);
  new_all_conds(3,:) = all_conds(3,:);
  new_all_conds(4,:) = all_conds_words(1,:);
  new_all_conds(5,:) = all_conds(6,:);
  
  all_conds = new_all_conds;
end
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
 
 
        
%      for x = 1:num_rest_TR_elim
%          all_conds(:,rest_ID_idx_run1(x)) = 99;
%          all_conds(:,rest_ID_idx_run2(x)) = 99;
%      end
%      
%      for r = 1:length(pre_conds)
%          if pre_conds(1,r) ~= 99;
%              trial_conds = pre_conds(:,r);
%              trial_runs = pre_runs(:,r);
%              trial_trials = pre_trials(:,r);
%              trial_TRs = pre_TRs(:,r);
%              all_conds = horzcat(all_conds, trial_conds);
%              all_runs = horzcat(all_runs, trial_runs);
%              all_trials = horzcat(all_trials,trial_trials);
%              all_TRs = horzcat(all_TRs, trial_TRs);
%          end
%          
%      end
%  end
 

  
  
      
%       9:424
%       my_cond = loc_conds(i);
%       my_run = loc_runs(i);
%       my_trial = loc_trial(i);
%       my_length = 26; 
%       %%%if first two subjects
%       %my_length = 18;
%       my_train_trs = 1:18;
%       %%%if first two subjects POSSIBLY
%       %my_train_trs = 9:18;
%       trial_runs = my_run*ones(1,my_length);
%       trial_nums = (my_trial+((my_run-1)*144))*ones(1,my_length);
%       trial_TR = 1:26;
%       trial_conds = zeros(count(conds_to_use),my_length);
%       if conds_to_use(my_cond)
%           trial_conds(my_cond, my_train_trs) = 1;
%       end
%       
%       
%       %if conds_to_use(my_cond)
%       all_conds = [all_conds trial_conds];
%       all_runs = [all_runs trial_runs];
%       all_trials = [all_trials trial_nums];
%       all_TRs = [all_TRs trial_TR];
%   end
% %   
% %   % prepare to insert 6s wait in all_runs
%   temp_run1 = ones(1,6);
%   temp_run2 = repmat(2,1,6);
%   temp_run3 = repmat(3,1,6);
  
  %prepare to insert 6s wait in all_conds
%   temp_conds = zeros(4,6);
%   
%   %prepare to insert 6s wait in all_trials 
%   temp_NaNs = NaN(1,6);
%   
%   %prepare to insert 6s wait in all_TRs
%   temp_TRs = 1:6;
%   % now i need to insert 6 seconds at the start of each run - temp matrix
%   % for conditions
%   all_conds = [temp_conds all_conds];
%   temp2_conds = all_conds(1:4,1:422); %Run 1 and includes 6s of wait. 
%   temp_endconds = all_conds(1:4,423:end); %Run 2 and 3 and does not include 6s of wait.
%   temp2_conds = [temp2_conds temp_conds]; %Run 1 and includes wait before and Run2 wait
%   all_conds = [temp2_conds temp_endconds]; %Run 1 and Run 2 with waits before and no wait between 2 and 3
%   temp2_conds = all_conds(1:4, 1:422);
%   temp3_conds = all_conds(1:4,423:844); %Run 2 - no waits
%   temp3_conds = [temp3_conds temp_conds]; %Run 2 with 6s wait after
%   temp_endconds = all_conds(1:4,845:end); %Run 3
%   all_conds = [temp2_conds temp3_conds temp_endconds]; 
%   
%   %for all_runs
%   all_runs = [temp_run1 all_runs];
%   temp2_runs = all_runs(1,1:422);
%   temp_endruns = all_runs(1,423:end);
%   temp2_runs = [temp2_runs temp_run2];
%   all_runs = [temp2_runs temp_endruns]; %Run 1 and Run 2 with waits before and no wait between 2 and 3
%   temp2_runs = all_runs(1, 1:422);
%   temp3_runs = all_runs(1,423:844); %Run 2 - no waits
%   temp3_runs = [temp3_runs temp_run3]; %Run 2 with 6s wait after
%   temp_endruns = all_runs(1,845:end); %Run 3
%   all_runs = [temp2_runs temp3_runs temp_endruns]; 
%   
%   %for all_trials
%   all_trials = [temp_NaNs all_trials];
%   temp2_trials = all_trials(1,1:422);
%   temp_endtrials = all_trials(1,423:end);
%   temp2_trials = [temp2_trials temp_NaNs];
%   all_trials = [temp2_trials temp_endtrials]; %Run 1 and Run 2 with waits before and no wait between 2 and 3
%   temp2_trials = all_trials(1, 1:422);
%   temp3_trials = all_trials(1,423:844); %Run 2 - no waits
%   temp3_trials = [temp3_trials temp_NaNs]; %Run 2 with 6s wait after
%   temp_endtrials = all_trials(1,845:end); %Run 3
%   all_trials = [temp2_trials temp3_trials temp_endtrials]; 
%   
%   
%   %for all_TRs
%   all_TRs = [temp_TRs all_TRs];
%   temp2_TRs = all_TRs(1,1:422);
%   temp_endTRs = all_TRs(1,423:end);
%   temp2_TRs = [temp2_TRs temp_NaNs];
%   all_TRs = [temp2_TRs temp_endTRs]; %Run 1 and Run 2 with waits before and no wait between 2 and 3
%   temp2_TRs = all_TRs(1, 1:422);
%   temp3_TRs = all_TRs(1,423:844); %Run 2 - no waits
%   temp3_TRs = [temp3_TRs temp_NaNs]; %Run 2 with 6s wait after
%   temp_endTRs = all_TRs(1,845:end); %Run 3
%   all_TRs = [temp2_TRs temp3_TRs temp_endTRs]; 
%   
  
% %   for i = 1:length(loc_conds)
% %     % get regressors for current trial
% %     my_cond = loc_conds(i);
% %     my_run = loc_runs(i);
% %     my_trial = loc_trial(i);
% %     my_length = 2; % 1.5s display stim and .5s iti    
% %     my_break = 8;  % break between miniblocks 
% %     my_startwait = 6;
% %     
% %     trial_runs = my_run*ones(1,my_length);
% %     trial_nums = (my_trial+((my_run-1)*144))*ones(1,my_length);
% %     trial_TR = 1:2;
% %     my_train_trs = 1:2; %training on both TRs of each trial.
% %     
%     create conditions matrix - build out regressors
%     trial_conds = zeros(length(conds_to_use),my_length);
% 
%      trial_conds = zeros(count(conds_to_use),my_length);
%     if conds_to_use(my_cond)
%         trial_conds(my_cond, my_train_trs) = 1;
%     end
    
% % %     % add in iti
% % %     if conds_to_use(end)
% % %       trial_conds(end, my_iti_trs) = 1;
% % %     end
% % %     
% % %     %if conds_to_use(my_cond)
% % %       all_conds = [all_conds trial_conds];
% % %       all_runs = [all_runs trial_runs];
% % %       all_trials = [all_trials trial_nums];
% % %       all_TRs = [all_TRs trial_TR];
% % %   end
% % %     

  
  fprintf('\n\n## using stimulus categories: %s\n',args.categories);
  whos all_regressors
  
  % add regressors & selectors to structure
  conds_name = sprintf('%s_conds',args.train_phase);
  runs_name = sprintf('%s_runs',args.train_phase);
  trials_name = sprintf('%s_trials',args.train_phase);
  TRs_name = sprintf('%s_TRs',args.train_phase);
  
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
  fprintf('\n\n## # of regressors in each condition: %s\n\n',mat2str(nansum(shifted_regs,2)'));
 
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
  
  fprintf('\n\n## reading in %s data from %s\n',args.train_phase, args.maskName);
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
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.train_phase, args.maskName, count(readmask));
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
  
  
  %%%TESTING PHASE
 

  %EmoDiF study - * 
  %**st_localizer - 2 blocks -213 TRs each block
  %**DF_encoding - 2 blocks -213 TRs each block
  %**Preview - 2 blocks - 183 TRs each block
  
  %mvpa_regs - %**1 = negative, 0 = neutral; 1= remember, 0 = forget
  %%%% READING IN REGRESSORS
  if strcmp(test_phase,'DFencode') == 1
      
  
  all_DF_runs = mvpa_regs.DFEncode.run;
  all_DF_trial = mvpa_regs.DFEncode.trial;
  all_DF_cat = mvpa_regs.DFEncode.cat;
  all_DF_TR = mvpa_regs.DFEncode.TR;
  all_DF_instr = mvpa_regs.DFEncode.instr;
  all_DF_emo = mvpa_regs.DFEncode.emo;
%   all_DF_subresp = mvpa_regs.DFEncode.subresp;
  
  
  
  %sets trial regressors by word category - Does it matter? maybe not.
  trial_regressors = zeros(length(args.categories),length(mvpa_regs.DFEncode.trial));
  
  %this stays empty.
  
  all_DF_regressors = trial_regressors;
    
  test_regs_name = sprintf('%s_all_regs',args.test_phase);
  test_runs_name = sprintf('%s_all_runs',args.test_phase);
  
  subj = initset_object(subj,'regressors',test_regs_name,all_DF_regressors);
  subj = initset_object(subj,'selector',test_runs_name,all_DF_runs);
  
   %%SHIFT REGRESSORS AND EXCLUDE REST TIMEPOINTS?
 
 
  % read in data from test phase masked by train phase voxels
  %
  fprintf('(+) study data for decoding\n');
  fprintf('(+) using study mask %s\n',args.maskName);
  
  %LOAD in test patterns
  
  cd (args.bold_dir);
  DF_raw_filenames = {'DF_encoding_1_corr_dt_mcf_brain.nii', 'DF_encoding_2_corr_dt_mcf_brain.nii'};
  
  readmask = get_mat(subj,'mask','read_mask');
  
  % read in EPI data as 'doubles'
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.test_phase, args.maskName, count(readmask));
  subj = load_spm_pattern(subj, 'DF_epis', 'read_mask', DF_raw_filenames);
  type ='double';
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.test_phase, args.maskName, count(readmask));
  
  
    %----------------------------------------------------------------------
  % z-score data
  %
  % Be sure to use the full 'runs' selector here so we include the 'rest'
  % TRs in the zscoring
  %
  fprintf('[+] z-scoring data\n\n');
  subj = zscore_runs(subj,'DF_epis',test_runs_name);
  
   %----------------------------------------------------------------------
  % setup decoding analysis!
  %
  train_epis_name = 'localizer_epis_z';
  test_epis_name = 'DF_epis_z';
  combined_epis_name = test_epis_name;
  mask_name = 'read_mask';
  
    % combine train & test epis
  decode_epis = [...
    get_masked_pattern(subj,train_epis_name,mask_name) ...
    get_mat(subj,'pattern',test_epis_name)];

  % don't need the train epis along anymore
  subj = remove_mat(subj,'pattern',train_epis_name);
  
  % replace the test_epis structure with the combined epis. this is important
  % so that the pattern has information of what mask it was created by.
  subj = set_mat(subj,'pattern',test_epis_name,decode_epis);
  
  % combine train & test phase regressors
  train_regs_name = sh_conds_name;
  combined_regs_name = sprintf('%s+%s',train_regs_name,test_regs_name);
  
  decode_regs = [ ...
    get_mat(subj,'regressors',train_regs_name) ...
    get_mat(subj,'regressors',test_regs_name)];
  
  subj = initset_object(subj,'regressors',combined_regs_name,decode_regs);
  
  % combine train & test phase selectors
  train_runs_name = norest_sels_name;
  combined_runs_name = sprintf('%s+%s',train_runs_name,test_runs_name);
  
  % set all testing runs to '2' for classifier decoding
  test_runs = get_mat(subj,'selector',test_runs_name);
  test_runs_for_decoding = 2*(test_runs>0);
  
  decode_runs = [...
    get_mat(subj,'selector',train_runs_name) ...
    test_runs_for_decoding];

  subj = initset_object(subj,'selector',combined_runs_name,decode_runs,...
    'group_name',combined_runs_name);  

  %-----------------------------------------------------------------------%
  % do decoding on combined data
  %
  
  switch(args.classifier)
      
      case 'bp'
          class_args.train_funct_name = 'train_bp';
          class_args.test_funct_name = 'test_bp';
          class_args.nHidden = 0;
          nIterations = args.repetitions;
          
      case 'L2logreg'
          
          class_args.train_funct_name = 'train_L2_RLR';
          class_args.test_funct_name = 'test_L2_RLR';
          class_args.penalty = args.penalty;
%           class_args.lambda = 'crossvalidation';
%           class_args.lambda = args.penalty;
          nIterations = 1;
          
      case 'logreg'
          class_args.train_funct_name = 'train_logreg';
          class_args.test_funct_name = 'test_logreg';
          class_args.penalty = str2double(penalty);
          
      case 'ridge'
          class_args.train_funct_name = 'train_ridge';
          class_args.test_funct_name = 'test_ridge';
          class_args.penalty = str2double(penalty);
          nIterations = 1;
          
      otherwise
          disp('(-) unknown classifier type!');
          return
          
  end
  
  for i = 1:nIterations
      fprintf('(+) decoding %s using %s (iter %d of %d)\n',...
          args.test_phase,args.train_phase,i,nIterations);
      
      [subj results] = ...
          cross_validation(subj, ...
          combined_epis_name, ...
          combined_regs_name, ...
          combined_runs_name, ...
          mask_name,...
          class_args);
      
      % save results structure for this iteration
     
      results_file = sprintf('%s/%s_%s_results.mat',...
           args.output_dir, args.subjID, args.test_phase);
      save(results_file, 'results')
      %write output to xls

      

      
      
      trial_info_test_perf=[all_DF_cat; all_DF_emo; all_DF_instr; all_DF_runs; all_DF_TR; all_DF_trial];
      
      dlmwrite(sprintf('%s/%s_%s_class_perf.txt', args.output_dir, args.subjID, args.test_phase),results.iterations.acts);
      dlmwrite(sprintf('%s/%s_%s.txt', args.output_dir, args.subjID, args.test_phase),trial_info_test_perf);
      
      fn = sprintf( '%s/%s_%s_parameters.txt', args.output_dir, args.subjID, args.test_phase);
      fid=fopen(fn,'w');
      fprintf(fid,sprintf('%s %s %s %s %s %s %s %s %s %s',args.subjID, args.train_phase, args.test_phase, maskName,  classifier, categories, penalty, shiftTRs,test_phase,rest_shift,word_combine));
      %     [rows cols] = size (concat_perf);
      %     x= repmat('%.4f\t',1,(cols-1));
      %     fprintf(fid,[x,'%.4f\n'],concat_perf);
      
      fclose(fid);
  end
  clear results 
  
  
  elseif strcmp(test_phase,'Preview') == 1
 
      all_preview_runs = mvpa_regs.preview.run;
      all_preview_cat = mvpa_regs.preview.cat;
      all_preview_trial = mvpa_regs.preview.trial;
      all_preview_instr = mvpa_regs.preview.instr;
%       all_preview_subresp = mvpa_regs.preview.subresp;
%       all_preview_rt = mvpa_regs.preview.rt;
      all_preview_TR = mvpa_regs.preview.TR;
      all_preview_emo = mvpa_regs.preview.emo;
      %put in preview regressors here. 
      
      
  
  
  %sets trial regressors by word category - Does it matter? maybe not.
  %sets trial regressors by word category - Does it matter? maybe not.
  trial_regressors = zeros(length(args.categories),length(mvpa_regs.DFEncode.trial));
  %this stays empty
  
  all_preview_regressors = trial_regressors;
    
  test_regs_name = sprintf('%s_all_regs',args.test_phase);
  test_runs_name = sprintf('%s_all_st_runs',args.test_phase);
  
  subj = initset_object(subj,'regressors',test_regs_name,all_preview_regressors);
  subj = initset_object(subj,'selector',test_runs_name,all_preview_runs);
  
   %%SHIFT REGRESSORS AND EXCLUDE REST TIMEPOINTS?
 
 
  % read in data from test phase masked by train phase voxels
  %
  fprintf('(+) study data for decoding\n');
  fprintf('(+) using study mask %s\n',args.maskName);
  
  %LOAD in test patterns
  
  cd (args.bold_dir);
  preview_raw_filenames = {'Preview1_corr_dt_mcf_brain.nii', 'Preview2_corr_dt_mcf_brain.nii'};
  
  readmask = get_mat(subj,'mask','read_mask');
  
  % read in EPI data as 'doubles'
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.test_phase, args.maskName, count(readmask));
  subj = load_spm_pattern(subj, 'preview_epis', 'read_mask', preview_raw_filenames);
  type ='double';
  fprintf('\n\n## loading %s data using %s with %d voxels\n', args.test_phase, args.maskName, count(readmask));
  
  
    %----------------------------------------------------------------------
  % z-score data
  %
  % Be sure to use the full 'runs' selector here so we include the 'rest'
  % TRs in the zscoring
  %
  fprintf('[+] z-scoring data\n\n');
  subj = zscore_runs(subj,'preview_epis',test_runs_name);
  
   %----------------------------------------------------------------------
  % setup decoding analysis!
  %
  train_epis_name = 'localizer_epis_z';
  test_epis_name = 'preview_epis_z';
  combined_epis_name = test_epis_name;
  mask_name = 'read_mask';
  
    % combine train & test epis
  decode_epis = [...
    get_masked_pattern(subj,train_epis_name,mask_name) ...
    get_mat(subj,'pattern',test_epis_name)];

  % don't need the train epis along anymore
  subj = remove_mat(subj,'pattern',train_epis_name);
  
  % replace the test_epis structure with the combined epis. this is important
  % so that the pattern has information of what mask it was created by.
  subj = set_mat(subj,'pattern',test_epis_name,decode_epis);
  
  % combine train & test phase regressors
  train_regs_name = sh_conds_name;
  combined_regs_name = sprintf('%s+%s',train_regs_name,test_regs_name);
  
  decode_regs = [ ...
    get_mat(subj,'regressors',train_regs_name) ...
    get_mat(subj,'regressors',test_regs_name)];
  
  subj = initset_object(subj,'regressors',combined_regs_name,decode_regs);
  
  % combine train & test phase selectors
  train_runs_name = norest_sels_name;
  combined_runs_name = sprintf('%s+%s',train_runs_name,test_runs_name);
  
  % set all testing runs to '2' for classifier decoding
  test_runs = get_mat(subj,'selector',test_runs_name);
  test_runs_for_decoding = 2*(test_runs>0);
  
  decode_runs = [...
    get_mat(subj,'selector',train_runs_name) ...
    test_runs_for_decoding];

  subj = initset_object(subj,'selector',combined_runs_name,decode_runs,...
    'group_name',combined_runs_name);  

  %-----------------------------------------------------------------------%
  % do decoding on combined data
  %
  
  switch(args.classifier)
      
      case 'bp'
          class_args.train_funct_name = 'train_bp';
          class_args.test_funct_name = 'test_bp';
          class_args.nHidden = 0;
          nIterations = args.repetitions;
          
      case 'L2logreg'
          
          class_args.train_funct_name = 'train_L2_RLR';
          class_args.test_funct_name = 'test_L2_RLR';
          class_args.penalty = args.penalty;
%           class_args.lambda = 'crossvalidation';
%           class_args.lambda = args.penalty;
          nIterations = 1;
          
      case 'logreg'
          class_args.train_funct_name = 'train_logreg';
          class_args.test_funct_name = 'test_logreg';
          class_args.penalty = str2double(penalty);
          
      case 'ridge'
          class_args.train_funct_name = 'train_ridge';
          class_args.test_funct_name = 'test_ridge';
          class_args.penalty = str2double(penalty);
          nIterations = 1;
          
      otherwise
          disp('(-) unknown classifier type!');
          return
          
  end
  
  for i = 1:nIterations
      fprintf('(+) decoding %s using %s (iter %d of %d)\n',...
          args.test_phase,args.train_phase,i,nIterations);
      
      [subj results] = ...
          cross_validation(subj, ...
          combined_epis_name, ...
          combined_regs_name, ...
          combined_runs_name, ...
          mask_name,...
          class_args);
      
      % save results structure for this iteration
     
      results_file = sprintf('%s/%s_%s_results.mat',...
           args.output_dir, args.subjID, args.test_phase);
      save(results_file, 'results')
      %write output to xls

      trial_info_test_perf=[all_preview_cat; all_preview_emo; all_preview_instr; all_preview_runs; all_preview_TR; all_preview_trial];
      
      dlmwrite(sprintf('%s/%s_%s_class_perf.txt', args.output_dir, args.subjID, args.test_phase),results.iterations.acts);
      dlmwrite(sprintf('%s/%s_%s_regressors.txt', args.output_dir, args.subjID, args.test_phase),trial_info_test_perf);
      
      fn = sprintf( '%s/%s_%s_parameters.txt', args.output_dir, args.subjID, args.test_phase);
      fid=fopen(fn,'w');
      fprintf(fid,sprintf('%s %s %s %s %s %s %s %s %s %s',args.subjID, args.train_phase, args.test_phase, maskName,  classifier, categories, penalty, shiftTRs,test_phase,rest_shift,word_combine));
      %     [rows cols] = size (concat_perf);
      %     x= repmat('%.4f\t',1,(cols-1));
      %     fprintf(fid,[x,'%.4f\n'],concat_perf);
      
      fclose(fid);
  end
  clear results 
  else
      
   disp('ERROR - no test_phase concurrance');
   
  end
      diary off
      cd(start_dir);
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
