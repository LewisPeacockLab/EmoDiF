function emodif_mvpa_decoding_parse(subjNum,test_phase,test_date)
%*emodif_mvpa_decoding_parse('101','DFencode','30-Apr-2018')
%test_phase can be 'DFencode' OR 'preview'
%SM can be 'hiconf', 'merged'
%doesn't need conditions or mask name 

  version = '2018Apr30';
  args.break = 3; %how many TRs are at the end of each run regardless of condition
  args.local.nTRs = 426; % 135 trials per run, 5 TRs per rest, each post miniblock. each single run 145 stim TRs and 75 rest TRs.  
  args.local.trialTR = 9;
  args.DF.nTRs = 426;
  args.DF.trialTR = 7; % word(1TR), blank(1TR), DFinstr(4TR), blank(1TR)
  args.DF.trialn = 30;
  args.DF.runs = 2;
  args.preview.nTRs = 366;
  args.preview.trialTR = 6; %word(1.5TR), scenes(1.5TRs), fixation(3TRs)
  args.preview.trialn = 30;
  args.preview.runs = 2;
  args.subjNum = subjNum;
  args.subjID = sprintf('emodif_%s',num2str(subjNum));
  
  
    %for astoria
   args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% % for tigger
%    args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);

  args.bold_dir = sprintf('%s/BOLD', args.subj_dir);
  args.mask_dir = sprintf('%s/mask', args.subj_dir);
  args.regs_dir = sprintf('%s/behav', args.subj_dir);
  args.output_dir = sprintf('%s/results/%s/%s',args.subj_dir, test_phase, test_date);

cd(args.output_dir)
if strcmp(test_phase, 'DFencode') == true
    class_perf = load(sprintf('emodif_%s_DFencode_class_perf.txt',subjNum));
    regressors = load(sprintf('emodif_%s_DFencode_regressors.txt',subjNum));
    cat = regressors(1,:);
    emo = regressors(2,:); %1 negative, 0 neutral
    instr = regressors(3,:); %1 remember, 0 forget
    runs = regressors(4,:);
    subresp = regressors(5,:);
    TR = regressors(6,:);
    trial = regressors(7,:);

elseif strcmp(test_phase, 'preview') == true
    class_perf = load(sprintf('emodif_%s_preview_class_perf.txt',subjNum));
    regressors = load(sprintf('emodif_%s_preview_regressors.txt',subjNum));
    cat = regressors(1,:);
    % RT = regressors(2,:);
    emo = regressors(2,:);
    instr = regressors(3,:);
    runs = regressors(4,:);
    subresp = regressors(5,:);
    TR = regressors(6,:);
    trial = regressors(7,:);
end

face_acts = class_perf(1,:);
scene_acts = class_perf(2,:);
object_acts = class_perf(3,:);
word_acts = class_perf(4,:);
rest_acts = class_perf(5,:);

all_regressors = vertcat(class_perf,regressors);

%start breaking everything up by trial 

if strcmp(test_phase, 'DFencode') == true
    instr_remember_idx = find(instr == 1);
    instr_forget_idx = find(instr == 0);
    instr_remember_neu_idx = find(instr == 1 & emo == 0);
    instr_remember_neg_idx = find(instr == 1 & emo == 1);
    instr_forget_neu_idx = find(instr == 0 & emo == 0);
    instr_forget_neg_idx = find(instr == 0 & emo == 1);
    
    %subsequent memory idx
    
    instr_forget_neu_1_idx = find(instr == 0 & emo == 0 & subresp == 1);
    instr_forget_neu_2_idx = find(instr == 0 & emo == 0 & subresp == 2);
    instr_forget_neu_3_idx = find(instr == 0 & emo == 0 & subresp == 3);
    instr_forget_neu_4_idx = find(instr == 0 & emo == 0 & subresp == 4);
    
    instr_forget_neg_1_idx = find(instr == 0 & emo == 1 & subresp == 1);
    instr_forget_neg_2_idx = find(instr == 0 & emo == 1 & subresp == 2);
    instr_forget_neg_3_idx = find(instr == 0 & emo == 1 & subresp == 3);    
    instr_forget_neg_4_idx = find(instr == 0 & emo == 1 & subresp == 4);
    
    instr_rem_neu_1_idx = find(instr == 1 & emo == 0 & subresp == 1);
    instr_rem_neu_2_idx = find(instr == 1 & emo == 0 & subresp == 2);
    instr_rem_neu_3_idx = find(instr == 1 & emo == 0 & subresp == 3);
    instr_rem_neu_4_idx = find(instr == 1 & emo == 0 & subresp == 4);
    
    instr_rem_neg_1_idx = find(instr == 1 & emo == 1 & subresp == 1);
    instr_rem_neg_2_idx = find(instr == 1 & emo == 1 & subresp == 2);
    instr_rem_neg_3_idx = find(instr == 1 & emo == 1 & subresp == 3);    
    instr_rem_neg_4_idx = find(instr == 1 & emo == 1 & subresp == 4);
    
    
    %check what type of subsequent memory
    
        %highconfidence
        instr_forget_neu_hiconf_OLD_idx = instr_forget_neu_1_idx;
        instr_forget_neg_hiconf_OLD_idx = instr_forget_neg_1_idx;
        instr_forget_hiconf_OLD_idx = horzcat(instr_forget_neu_hiconf_OLD_idx, instr_forget_neg_hiconf_OLD_idx);
        
        instr_forget_neu_hiconf_NEW_idx = horzcat(instr_forget_neu_2_idx, instr_forget_neu_3_idx, instr_forget_neu_4_idx);
        instr_forget_neg_hiconf_NEW_idx = horzcat(instr_forget_neg_2_idx, instr_forget_neg_3_idx, instr_forget_neg_4_idx);
        instr_forget_hiconf_NEW_idx = horzcat(instr_forget_neu_hiconf_NEW_idx, instr_forget_neg_hiconf_NEW_idx);
        
        instr_rem_neu_hiconf_OLD_idx = instr_rem_neu_1_idx;
        instr_rem_neg_hiconf_OLD_idx = instr_rem_neg_1_idx;
        instr_rem_hiconf_OLD_idx = horzcat(instr_rem_neu_hiconf_OLD_idx, instr_rem_neg_hiconf_OLD_idx);
        
        instr_rem_neu_hiconf_NEW_idx = horzcat(instr_rem_neu_2_idx, instr_rem_neu_3_idx, instr_rem_neu_4_idx);
        instr_rem_neg_hiconf_NEW_idx = horzcat(instr_rem_neg_2_idx, instr_rem_neg_3_idx, instr_rem_neg_4_idx);
        instr_rem_hiconf_NEW_idx = horzcat(instr_rem_neu_hiconf_NEW_idx, instr_rem_neg_hiconf_NEW_idx);
        
        
        
        % high confidence means only high confidence old is counted as old,
        % all other responses count as new
        %merged
        instr_forget_neu_merged_OLD_idx = horzcat(instr_forget_neu_1_idx, instr_forget_neu_2_idx);
        instr_forget_neg_merged_OLD_idx = horzcat(instr_forget_neg_1_idx, instr_forget_neg_2_idx);
        instr_forget_merged_OLD_idx = horzcat(instr_forget_neu_merged_OLD_idx, instr_forget_neg_merged_OLD_idx);
        
        instr_forget_neu_merged_NEW_idx = horzcat(instr_forget_neu_3_idx, instr_forget_neu_4_idx);
        instr_forget_neg_merged_NEW_idx = horzcat(instr_forget_neg_3_idx, instr_forget_neg_4_idx);
        instr_forget_merged_NEW_idx = horzcat(instr_forget_neu_merged_NEW_idx, instr_forget_neg_merged_NEW_idx);
        
        instr_rem_neu_merged_OLD_idx = horzcat(instr_rem_neu_1_idx, instr_rem_neu_2_idx);
        instr_rem_neg_merged_OLD_idx = horzcat(instr_rem_neg_1_idx, instr_rem_neg_2_idx);
        instr_rem_merged_OLD_idx = horzcat(instr_rem_neu_merged_OLD_idx, instr_rem_neg_merged_OLD_idx);
        
        instr_rem_neu_merged_NEW_idx = horzcat(instr_rem_neu_3_idx, instr_rem_neu_4_idx);
        instr_rem_neg_merged_NEW_idx = horzcat(instr_rem_neg_3_idx, instr_rem_neg_4_idx);
        instr_rem_merged_NEW_idx = horzcat(instr_rem_neu_merged_NEW_idx, instr_rem_neg_merged_NEW_idx);

         % merged means that confidence levels are merged so that
         % confident old + unconfident old = old
         % confident new + unconfident new = new
    
    

%     instr_remember_neu_idx = find(instr_remember_idx
    
%     for k = 1:length(args.DF.trialn)
%         remember_reg.face(1,(k-1):() 

    %this is where we should put the rest of the regressors as a sanity
    %cut regressors for 1 ST 
    %check
    
    for i = 1:length(instr_remember_idx)
            remember_acts(:,i)=class_perf(:,instr_remember_idx(i));
    end
        
rem_face = reshape(remember_acts(1,:),7,30)';
rem.all.acts.face = rem_face;
rem_scene =reshape(remember_acts(2,:),7,30)';
rem.all.acts.scene = rem_scene;
rem_object =reshape(remember_acts(3,:),7,30)';
rem.all.acts.object = rem_object;
rem_word =reshape(remember_acts(4,:),7,30)';
rem.all.acts.word = rem_word;
rem_rest =   reshape(remember_acts(5,:),7,30)';
rem.all.acts.rest = rem_rest;
rem_check = reshape(instr_remember_idx,7,30)';
rem.all.acts.check = rem_check;


    for i = 1:length(instr_forget_idx)
            forget_acts(:,i)=class_perf(:,instr_forget_idx(i));
    end
        
for_face = reshape(forget_acts(1,:),7,30)';
forget.all.acts.face = for_face;
for_scene =reshape(forget_acts(2,:),7,30)';
forget.all.acts.scene = for_scene;
for_object =reshape(forget_acts(3,:),7,30)';
forget.all.acts.object = for_object;
for_word =reshape(forget_acts(4,:),7,30)';
forget.all.acts.word = for_word;
for_rest =   reshape(forget_acts(5,:),7,30)';
forget.all.acts.rest = for_rest;
for_check = reshape(instr_remember_idx,7,30)';
forget.all.acts.check = for_check;

    for i = 1:length(instr_remember_neu_idx)
            remember_neu_acts(:,i)=class_perf(:,instr_remember_neu_idx(i));
    end

rem_neu_face = reshape(remember_neu_acts(1,:),7,15)';
rem.neu.acts.face = rem_neu_face;
rem_neu_scene =reshape(remember_neu_acts(2,:),7,15)';
rem.neu.acts.scene = rem_neu_scene;
rem_neu_object =reshape(remember_neu_acts(3,:),7,15)';
rem.neu.acts.object = rem_neu_object;
rem_neu_word =reshape(remember_neu_acts(4,:),7,15)';
rem.neu.acts.word = rem_neu_word;
rem_neu_rest =   reshape(remember_neu_acts(5,:),7,15)';
rem.neu.acts.rest = rem_neu_rest;
rem_neu_check = reshape(instr_remember_neu_idx,7,15)';
rem.neu.acts.check = rem_neu_check;

    for i = 1:length(instr_remember_neg_idx)
            remember_neg_acts(:,i)=class_perf(:,instr_remember_neg_idx(i));
    end

rem_neg_face = reshape(remember_neg_acts(1,:),7,15)';
rem.neg.acts.face = rem_neg_face;
rem_neg_scene =reshape(remember_neg_acts(2,:),7,15)';
rem.neg.acts.scene = rem_neg_scene;
rem_neg_object =reshape(remember_neg_acts(3,:),7,15)';
rem.neg.acts.object = rem_neg_object;
rem_neg_word =reshape(remember_neg_acts(4,:),7,15)';
rem.neg.acts.word = rem_neg_word;
rem_neg_rest =   reshape(remember_neg_acts(5,:),7,15)';
rem.neg.acts.rest = rem_neg_rest;
rem_neg_check = reshape(instr_remember_neg_idx,7,15)';
rem.neg.acts.check = rem_neg_check;

    for i = 1:length(instr_forget_neu_idx)
            forget_neu_acts(:,i)=class_perf(:,instr_forget_neu_idx(i));
    end

forget_neu_face = reshape(forget_neu_acts(1,:),7,15)';
forget.neu.acts.face = forget_neu_face;
forget_neu_scene =reshape(forget_neu_acts(2,:),7,15)';
forget.neu.acts.scene = forget_neu_scene;
forget_neu_object =reshape(forget_neu_acts(3,:),7,15)';
forget.neu.acts.object = forget_neu_object;
forget_neu_word =reshape(forget_neu_acts(4,:),7,15)';
forget.neu.acts.word = forget_neu_word;
forget_neu_rest =   reshape(forget_neu_acts(5,:),7,15)';
forget.neu.acts.rest = forget_neu_rest;
forget_neu_check = reshape(instr_forget_neu_idx,7,15)';
forget.neu.acts.check = forget_neu_check;

    for i = 1:length(instr_forget_neg_idx)
            forget_neg_acts(:,i)=class_perf(:,instr_forget_neg_idx(i));
    end

forget_neg_face = reshape(forget_neg_acts(1,:),7,15)';
forget.neg.acts.face = forget_neg_face;
forget_neg_scene =reshape(forget_neg_acts(2,:),7,15)';
forget.neg.acts.scene = forget_neg_scene;
forget_neg_object =reshape(forget_neg_acts(3,:),7,15)';
forget.neg.acts.object = forget_neg_object;
forget_neg_word =reshape(forget_neg_acts(4,:),7,15)';
forget.neg.acts.word = forget_neg_word;
forget_neg_rest =   reshape(forget_neg_acts(5,:),7,15)';
forget.neg.acts.rest = forget_neg_rest;
forget_neg_check = reshape(instr_forget_neg_idx,7,15)';
forget.neg.acts.check = forget_neg_check;

    for i = 1:length(instr_forget_hiconf_OLD_idx)
            forget_hiconf_OLD_acts(:,i)=class_perf(:,instr_forget_hiconf_OLD_idx(i));
    end
 
    %high confidence SM
forget_hiconf_OLD_n = length(instr_forget_hiconf_OLD_idx)/7;
forget.highconf.acts.OLD.num = forget_hiconf_OLD_n;
forget_hiconf_OLD_face = reshape(forget_hiconf_OLD_acts(1,:),7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.face = forget_hiconf_OLD_face;
forget_hiconf_OLD_scene = reshape(forget_hiconf_OLD_acts(2,:),7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.scene = forget_hiconf_OLD_scene;
forget_hiconf_OLD_object = reshape(forget_hiconf_OLD_acts(3,:),7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.object = forget_hiconf_OLD_object;
forget_hiconf_OLD_word = reshape(forget_hiconf_OLD_acts(4,:),7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.word = forget_hiconf_OLD_word;
forget_hiconf_OLD_rest = reshape(forget_hiconf_OLD_acts(5,:),7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.rest = forget_hiconf_OLD_rest;
instr_forget_hiconf_OLD_check = reshape(instr_forget_hiconf_OLD_idx,7,forget_hiconf_OLD_n)';
forget.highconf.acts.OLD.check = instr_forget_hiconf_OLD_check;
    
for i = 1:length(instr_forget_hiconf_NEW_idx)
    forget_hiconf_NEW_acts(:,i)=class_perf(:,instr_forget_hiconf_NEW_idx(i));
end

forget_hiconf_NEW_n = length(instr_forget_hiconf_NEW_idx)/7;
forget.highconf.acts.NEW.num = forget_hiconf_NEW_n;
forget_hiconf_NEW_face = reshape(forget_hiconf_NEW_acts(1,:),7,forget_hiconf_NEW_n)';
forget.highconf.acts.NEW.face = forget_hiconf_NEW_face;
forget_hiconf_NEW_scene = reshape(forget_hiconf_NEW_acts(2,:),7,forget_hiconf_NEW_n');
forget.highconf.acts.NEW.scene = forget_hiconf_NEW_scene;
forget_hiconf_NEW_object = reshape(forget_hiconf_NEW_acts(3,:),7,forget_hiconf_NEW_n)';
forget.highconf.acts.NEW.object = forget_hiconf_NEW_object;
forget_hiconf_NEW_word = reshape(forget_hiconf_NEW_acts(4,:),7,forget_hiconf_NEW_n)';
forget.highconf.acts.NEW.word = forget_hiconf_NEW_word;
forget_hiconf_NEW_rest = reshape(forget_hiconf_NEW_acts(5,:),7,forget_hiconf_NEW_n)';
forget.highconf.acts.NEW.rest = forget_hiconf_NEW_rest;
instr_forget_hiconf_NEW_check = reshape(instr_forget_hiconf_NEW_idx,7,forget_hiconf_NEW_n)';
forget.highconf.acts.NEW.check = instr_forget_hiconf_NEW_check;

for i = 1:length(instr_rem_hiconf_OLD_idx)
    rem_hiconf_OLD_acts(:,i)=class_perf(:,instr_rem_hiconf_OLD_idx(i));
end
    
rem_hiconf_OLD_n = length(instr_rem_hiconf_OLD_idx)/7;
rem.highconf.acts.OLD.num = rem_hiconf_OLD_n;
rem_hiconf_OLD_face = reshape(rem_hiconf_OLD_acts(1,:),7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.face = rem_hiconf_OLD_face;
rem_hiconf_OLD_scene = reshape(rem_hiconf_OLD_acts(2,:),7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.scene = rem_hiconf_OLD_scene;
rem_hiconf_OLD_object = reshape(rem_hiconf_OLD_acts(3,:),7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.object = rem_hiconf_OLD_object;
rem_hiconf_OLD_word = reshape(rem_hiconf_OLD_acts(4,:),7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.word = rem_hiconf_OLD_word;
rem_hiconf_OLD_rest = reshape(rem_hiconf_OLD_acts(5,:),7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.rest = rem_hiconf_OLD_rest;
instr_rem_hiconf_OLD_check = reshape(instr_rem_hiconf_OLD_idx,7,rem_hiconf_OLD_n)';
rem.highconf.acts.OLD.check = instr_rem_hiconf_OLD_check;

    
for i = 1:length(instr_rem_hiconf_NEW_idx)
    rem_hiconf_NEW_acts(:,i)=class_perf(:,instr_rem_hiconf_NEW_idx(i));
end

rem_hiconf_NEW_n = length(instr_rem_hiconf_NEW_idx)/7;
rem.highconf.acts.NEW.num= rem_hiconf_NEW_n;
rem_hiconf_NEW_face = reshape(rem_hiconf_NEW_acts(1,:),7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.face = rem_hiconf_NEW_face;
rem_hiconf_NEW_scene = reshape(rem_hiconf_NEW_acts(2,:),7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.scene = rem_hiconf_NEW_scene;
rem_hiconf_NEW_object = reshape(rem_hiconf_NEW_acts(3,:),7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.object = rem_hiconf_NEW_object;
rem_hiconf_NEW_word = reshape(rem_hiconf_NEW_acts(4,:),7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.word = rem_hiconf_NEW_word;
rem_hiconf_NEW_rest = reshape(rem_hiconf_NEW_acts(5,:),7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.rest = rem_hiconf_NEW_rest;
instr_rem_hiconf_NEW_check = reshape(instr_rem_hiconf_NEW_idx,7,rem_hiconf_NEW_n)';
rem.highconf.acts.NEW.check = instr_rem_hiconf_NEW_check;

  % MERGED confidence SM
  
for i = 1:length(instr_forget_merged_OLD_idx)
    forget_merged_OLD_acts(:,i)=class_perf(:,instr_forget_merged_OLD_idx(i));
end

forget_merged_OLD_n = length(instr_forget_merged_OLD_idx)/7;
forget.merged.acts.OLD.num = forget_merged_OLD_n;
forget_merged_OLD_face = reshape(forget_merged_OLD_acts(1,:),7,forget_merged_OLD_n)';
forget.merged.acts.OLD.face = forget_merged_OLD_face;
forget_merged_OLD_scene = reshape(forget_merged_OLD_acts(2,:),7,forget_merged_OLD_n)';
forget.merged.acts.OLD.scene = forget_merged_OLD_scene;
forget_merged_OLD_object = reshape(forget_merged_OLD_acts(3,:),7,forget_merged_OLD_n)';
forget.merged.acts.OLD.object = forget_merged_OLD_object;
forget_merged_OLD_word = reshape(forget_merged_OLD_acts(4,:),7,forget_merged_OLD_n)';
forget.merged.acts.OLD.word = forget_merged_OLD_word;
forget_merged_OLD_rest = reshape(forget_merged_OLD_acts(5,:),7,forget_merged_OLD_n)';
forget.merged.acts.OLD.rest = forget_merged_OLD_rest;
instr_forget_merged_OLD_check = reshape(instr_forget_merged_OLD_idx,7,forget_merged_OLD_n)';
forget.merged.acts.OLD.check = instr_forget_merged_OLD_check;

for i = 1:length(instr_forget_merged_NEW_idx)
    forget_merged_NEW_acts(:,i)=class_perf(:,instr_forget_merged_NEW_idx(i));
end

forget_merged_NEW_n = length(instr_forget_merged_NEW_idx)/7;
forget.merged.acts.NEW.num = forget_merged_NEW_n;
forget_merged_NEW_face = reshape(forget_merged_NEW_acts(1,:),7,forget_merged_NEW_n)';
forget.merged.acts.NEW.face = forget_merged_NEW_face;
forget_merged_NEW_scene = reshape(forget_merged_NEW_acts(2,:),7,forget_merged_NEW_n)';
forget.merged.acts.NEW.scene = forget_merged_NEW_scene;
forget_merged_NEW_object = reshape(forget_merged_NEW_acts(3,:),7,forget_merged_NEW_n)';
forget.merged.acts.NEW.object = forget_merged_NEW_object;
forget_merged_NEW_word = reshape(forget_merged_NEW_acts(4,:),7,forget_merged_NEW_n)';
forget.merged.acts.NEW.word = forget_merged_NEW_word;
forget_merged_NEW_rest = reshape(forget_merged_NEW_acts(5,:),7,forget_merged_NEW_n)';
forget.merged.acts.NEW.rest = forget_merged_NEW_rest;
instr_forget_merged_NEW_check = reshape(instr_forget_merged_NEW_idx,7,forget_merged_NEW_n)';
forget.merged.acts.NEW.check = instr_forget_merged_NEW_check;

for i = 1:length(instr_rem_merged_OLD_idx)
    rem_merged_OLD_acts(:,i)=class_perf(:,instr_rem_merged_OLD_idx(i));
end

rem_merged_OLD_n = length(instr_rem_merged_OLD_idx)/7;
rem.merged.acts.OLD.num = rem_merged_OLD_n;
rem_merged_OLD_face = reshape(rem_merged_OLD_acts(1,:),7,rem_merged_OLD_n)';
rem.merged.acts.OLD.face = rem_merged_OLD_face;
rem_merged_OLD_scene = reshape(rem_merged_OLD_acts(2,:),7,rem_merged_OLD_n)';
rem.merged.acts.OLD.scene = rem_merged_OLD_scene;
rem_merged_OLD_object = reshape(rem_merged_OLD_acts(3,:),7,rem_merged_OLD_n)';
rem.merged.acts.OLD.object = rem_merged_OLD_object;
rem_merged_OLD_word = reshape(rem_merged_OLD_acts(4,:),7,rem_merged_OLD_n)';
rem.merged.acts.OLD.word = rem_merged_OLD_word;
rem_merged_OLD_rest = reshape(rem_merged_OLD_acts(5,:),7,rem_merged_OLD_n)';
rem.merged.acts.OLD.rest = rem_merged_OLD_rest;
instr_rem_merged_OLD_check = reshape(instr_rem_merged_OLD_idx,7,rem_merged_OLD_n)';
rem.merged.acts.OLD.check = instr_rem_merged_OLD_check;


for i = 1:length(instr_rem_merged_NEW_idx)
    rem_merged_NEW_acts(:,i)=class_perf(:,instr_rem_merged_NEW_idx(i));
end

rem_merged_NEW_n = length(instr_rem_merged_NEW_idx)/7;
rem.merged.acts.NEW.num= rem_merged_NEW_n;
rem_merged_NEW_face = reshape(rem_merged_NEW_acts(1,:),7,rem_merged_NEW_n)';
rem.merged.acts.NEW.face = rem_merged_NEW_face;
rem_merged_NEW_scene = reshape(rem_merged_NEW_acts(2,:),7,rem_merged_NEW_n)';
rem.merged.acts.NEW.scene = rem_merged_NEW_scene;
rem_merged_NEW_object = reshape(rem_merged_NEW_acts(3,:),7,rem_merged_NEW_n)';
rem.merged.acts.NEW.object = rem_merged_NEW_object;
rem_merged_NEW_word = reshape(rem_merged_NEW_acts(4,:),7,rem_merged_NEW_n)';
rem.merged.acts.NEW.word = rem_merged_NEW_word;
rem_merged_NEW_rest = reshape(rem_merged_NEW_acts(5,:),7,rem_merged_NEW_n)';
rem.merged.acts.NEW.rest = rem_merged_NEW_rest;
instr_rem_merged_NEW_check = reshape(instr_rem_merged_NEW_idx,7,rem_merged_NEW_n)';
rem.merged.acts.NEW.check = instr_rem_merged_NEW_check;



%%% averaging %%%

rem.all.mean.face = mean(rem_face);
rem.all.mean.scene = mean(rem_scene);
rem.all.mean.object = mean(rem_object);
rem.all.mean.word = mean(rem_word);
rem.all.mean.rest = mean(rem_rest);

forget.all.mean.face = mean(for_face);
forget.all.mean.scene = mean(for_scene);
forget.all.mean.object = mean(for_object);
forget.all.mean.word = mean(for_word);
forget.all.mean.rest = mean(for_rest);

rem.neu.mean.face = mean(rem_neu_face);
rem.neu.mean.scene = mean(rem_neu_scene);
rem.neu.mean.object = mean(rem_neu_object);
rem.neu.mean.word = mean(rem_neu_word);
rem.neu.mean.rest = mean(rem_neu_rest);

rem.neg.mean.face = mean(rem_neg_face);
rem.neg.mean.scene = mean(rem_neg_scene);
rem.neg.mean.object = mean(rem_neg_object);
rem.neg.mean.word = mean(rem_neg_word);
rem.neg.mean.rest = mean(rem_neg_rest);

forget.neu.mean.face = mean(forget_neu_face);
forget.neu.mean.scene = mean(forget_neu_scene);
forget.neu.mean.object = mean(forget_neu_object);
forget.neu.mean.word = mean(forget_neu_word);
forget.neu.mean.rest = mean(forget_neu_rest);

forget.neg.mean.face = mean(forget_neg_face);
forget.neg.mean.scene = mean(forget_neg_scene);
forget.neg.mean.object = mean(forget_neg_object);
forget.neg.mean.word = mean(forget_neg_word);
forget.neg.mean.rest = mean(forget_neg_rest);

%now by SM

forget.highconf.mean.OLD.face = mean(forget_hiconf_OLD_face);
forget.highconf.mean.OLD.scene= mean(forget_hiconf_OLD_scene);
forget.highconf.mean.OLD.object= mean(forget_hiconf_OLD_object);
forget.highconf.mean.OLD.word= mean(forget_hiconf_OLD_word);
forget.highconf.mean.OLD.rest= mean(forget_hiconf_OLD_rest);

forget.highconf.mean.NEW.face = mean(forget_hiconf_NEW_face);
forget.highconf.mean.NEW.scene= mean(forget_hiconf_NEW_scene);
forget.highconf.mean.NEW.object= mean(forget_hiconf_NEW_object);
forget.highconf.mean.NEW.word= mean(forget_hiconf_NEW_word);
forget.highconf.mean.NEW.rest= mean(forget_hiconf_NEW_rest);

rem.highconf.mean.OLD.face = mean(rem_hiconf_OLD_face);
rem.highconf.mean.OLD.scene= mean(rem_hiconf_OLD_scene);
rem.highconf.mean.OLD.object= mean(rem_hiconf_OLD_object);
rem.highconf.mean.OLD.word= mean(rem_hiconf_OLD_word);
rem.highconf.mean.OLD.rest= mean(rem_hiconf_OLD_rest);

rem.highconf.mean.NEW.face = mean(rem_hiconf_NEW_face);
rem.highconf.mean.NEW.scene= mean(rem_hiconf_NEW_scene);
rem.highconf.mean.NEW.object= mean(rem_hiconf_NEW_object);
rem.highconf.mean.NEW.word= mean(rem_hiconf_NEW_word);
rem.highconf.mean.NEW.rest= mean(rem_hiconf_NEW_rest);

forget.merged.mean.OLD.face = mean(forget_hiconf_OLD_face);
forget.merged.mean.OLD.scene= mean(forget_hiconf_OLD_scene);
forget.merged.mean.OLD.object= mean(forget_hiconf_OLD_object);
forget.merged.mean.OLD.word= mean(forget_hiconf_OLD_word);
forget.merged.mean.OLD.rest= mean(forget_hiconf_OLD_rest);

forget.merged.mean.NEW.face = mean(forget_hiconf_NEW_face);
forget.merged.mean.NEW.scene= mean(forget_hiconf_NEW_scene);
forget.merged.mean.NEW.object= mean(forget_hiconf_NEW_object);
forget.merged.mean.NEW.word= mean(forget_hiconf_NEW_word);
forget.merged.mean.NEW.rest= mean(forget_hiconf_NEW_rest);

rem.merged.mean.OLD.face = mean(rem_hiconf_OLD_face);
rem.merged.mean.OLD.scene= mean(rem_hiconf_OLD_scene);
rem.merged.mean.OLD.object= mean(rem_hiconf_OLD_object);
rem.merged.mean.OLD.word= mean(rem_hiconf_OLD_word);
rem.merged.mean.OLD.rest= mean(rem_hiconf_OLD_rest);

rem.merged.mean.NEW.face = mean(rem_hiconf_NEW_face);
rem.merged.mean.NEW.scene= mean(rem_hiconf_NEW_scene);
rem.merged.mean.NEW.object= mean(rem_hiconf_NEW_object);
rem.merged.mean.NEW.word= mean(rem_hiconf_NEW_word);
rem.merged.mean.NEW.rest= mean(rem_hiconf_NEW_rest);




%%standard error%%

rem.all.std.face = std(rem_face);
rem.all.std.scene = std(rem_scene);
rem.all.std.object = std(rem_object);
rem.all.std.word = std(rem_word);
rem.all.std.rest = std(rem_rest);
   
forget.all.std.face = std(for_face);
forget.all.std.scene = std(for_scene);
forget.all.std.object = std(for_object);
forget.all.std.word = std(for_word);
forget.all.std.rest = std(for_rest);

rem.neu.std.face = std(rem_neu_face);
rem.neu.std.scene = std(rem_neu_scene);
rem.neu.std.object = std(rem_neu_object);
rem.neu.std.word = std(rem_neu_word);
rem.neu.std.rest = std(rem_neu_rest);

rem.neg.std.face = std(rem_neg_face);
rem.neg.std.scene = std(rem_neg_scene);
rem.neg.std.object = std(rem_neg_object);
rem.neg.std.word = std(rem_neg_word);
rem.neg.std.rest = std(rem_neg_rest);

forget.neu.std.face = std(forget_neu_face);
forget.neu.std.scene = std(forget_neu_scene);
forget.neu.std.object = std(forget_neu_object);
forget.neu.std.word = std(forget_neu_word);
forget.neu.std.rest = std(forget_neu_rest);

forget.neg.std.face = std(forget_neg_face);
forget.neg.std.scene = std(forget_neg_scene);
forget.neg.std.object = std(forget_neg_object);
forget.neg.std.word = std(forget_neg_word);
forget.neg.std.rest = std(forget_neg_rest);

%now by std

forget.highconf.std.OLD.face = std(forget_hiconf_OLD_face);
forget.highconf.std.OLD.scene= std(forget_hiconf_OLD_scene);
forget.highconf.std.OLD.object= std(forget_hiconf_OLD_object);
forget.highconf.std.OLD.word= std(forget_hiconf_OLD_word);
forget.highconf.std.OLD.rest= std(forget_hiconf_OLD_rest);

forget.highconf.std.NEW.face = std(forget_hiconf_NEW_face);
forget.highconf.std.NEW.scene= std(forget_hiconf_NEW_scene);
forget.highconf.std.NEW.object= std(forget_hiconf_NEW_object);
forget.highconf.std.NEW.word= std(forget_hiconf_NEW_word);
forget.highconf.std.NEW.rest= std(forget_hiconf_NEW_rest);

rem.highconf.std.OLD.face = std(rem_hiconf_OLD_face);
rem.highconf.std.OLD.scene= std(rem_hiconf_OLD_scene);
rem.highconf.std.OLD.object= std(rem_hiconf_OLD_object);
rem.highconf.std.OLD.word= std(rem_hiconf_OLD_word);
rem.highconf.std.OLD.rest= std(rem_hiconf_OLD_rest);

rem.highconf.std.NEW.face = std(rem_hiconf_NEW_face);
rem.highconf.std.NEW.scene= std(rem_hiconf_NEW_scene);
rem.highconf.std.NEW.object= std(rem_hiconf_NEW_object);
rem.highconf.std.NEW.word= std(rem_hiconf_NEW_word);
rem.highconf.std.NEW.rest= std(rem_hiconf_NEW_rest);

forget.merged.std.OLD.face = std(forget_hiconf_OLD_face);
forget.merged.std.OLD.scene= std(forget_hiconf_OLD_scene);
forget.merged.std.OLD.object= std(forget_hiconf_OLD_object);
forget.merged.std.OLD.word= std(forget_hiconf_OLD_word);
forget.merged.std.OLD.rest= std(forget_hiconf_OLD_rest);

forget.merged.std.NEW.face = std(forget_hiconf_NEW_face);
forget.merged.std.NEW.scene= std(forget_hiconf_NEW_scene);
forget.merged.std.NEW.object= std(forget_hiconf_NEW_object);
forget.merged.std.NEW.word= std(forget_hiconf_NEW_word);
forget.merged.std.NEW.rest= std(forget_hiconf_NEW_rest);

rem.merged.std.OLD.face = std(rem_hiconf_OLD_face);
rem.merged.std.OLD.scene= std(rem_hiconf_OLD_scene);
rem.merged.std.OLD.object= std(rem_hiconf_OLD_object);
rem.merged.std.OLD.word= std(rem_hiconf_OLD_word);
rem.merged.std.OLD.rest= std(rem_hiconf_OLD_rest);

rem.merged.std.NEW.face = std(rem_hiconf_NEW_face);
rem.merged.std.NEW.scene= std(rem_hiconf_NEW_scene);
rem.merged.std.NEW.object= std(rem_hiconf_NEW_object);
rem.merged.std.NEW.word= std(rem_hiconf_NEW_word);
rem.merged.std.NEW.rest= std(rem_hiconf_NEW_rest);

results.remember = rem;
results.forget = forget;



filename = sprintf('results_%s.mat',test_phase);
save(filename,'results');

%%%% back to preview %%%%%%
elseif strcmp(test_phase, 'preview') == true;
    emo_negative_idx = find(emo == 1);
    emo_neutral_idx = find(emo == 0);
    
    for i = 1:length(emo_negative_idx)
            negative_acts(:,i)=class_perf(:,emo_negative_idx(i));
    end
        
neg_face = reshape(negative_acts(1,:),6,30)';
neg_scene =reshape(negative_acts(2,:),6,30)';
neg_object =reshape(negative_acts(3,:),6,30)';
neg_word =reshape(negative_acts(4,:),6,30)';
neg_rest =   reshape(negative_acts(5,:),6,30)';

    for i = 1:length(emo_neutral_idx)
            neutral_acts(:,i)=class_perf(:,emo_neutral_idx(i));
    end
        
neu_face = reshape(neutral_acts(1,:),6,30)';
neutral.acts.face = neu_face;
neu_scene =reshape(neutral_acts(2,:),6,30)';
neutral.acts.scene = neu_scene;
neu_object =reshape(neutral_acts(3,:),6,30)';
neutral.acts.object = neu_object;
neu_word =reshape(neutral_acts(4,:),6,30)';
neutral.acts.word = neu_word;
neu_rest =   reshape(neutral_acts(5,:),6,30)';
neutral.acts.rest = neu_rest;
            
negative.mean.face = mean(neg_face);
negative.mean.scene = mean(neg_scene);
negative.mean.object = mean(neg_object);
negative.mean.word = mean(neg_word);
negative.mean.rest = mean(neg_rest);
   
neutral.mean.face = mean(neu_face);
neutral.mean.scene = mean(neu_scene);
neutral.mean.object = mean(neu_object);
neutral.mean.word = mean(neu_word);
neutral.mean.rest = mean(neu_rest);

negative.std.face = std(neg_face);
negative.std.scene = std(neg_scene);
negative.std.object = std(neg_object);
negative.std.word = std(neg_word);
negative.std.rest = std(neg_rest);
   
neutral.std.face = std(neu_face);
neutral.std.scene = std(neu_scene);
neutral.std.object = std(neu_object);
neutral.std.word = std(neu_word);
neutral.std.rest = std(neu_rest);

results.preview.negative = negative;
results.preview.neutral = neutral;

filename = sprintf('results_%s.mat',test_phase);
save(filename,'results');

end







end