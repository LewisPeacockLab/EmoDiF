function emodif_mvpa_decoding_parse(subjNum,maskName,conditions,test_phase, test_date);
%*emodif_mvpa_decoding_parse('101','tempoccfusi_pHg_LOC_combined_epi_space','fsoner', 'DFencode','30-Apr-2018')
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
%    args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
   args.subj_dir = sprintf('/Users/TWang/emodif_data/%s', args.subjID);

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
rem.acts.face = rem_face;
rem_scene =reshape(remember_acts(2,:),7,30)';
rem.acts.scene = rem_scene;
rem_object =reshape(remember_acts(3,:),7,30)';
rem.acts.object = rem_object;
rem_word =reshape(remember_acts(4,:),7,30)';
rem.acts.word = rem_word;
rem_rest =   reshape(remember_acts(5,:),7,30)';
rem.acts.rest = rem_rest;
rem_check = reshape(instr_remember_idx,7,30)';
rem.acts.check = rem_check;


    for i = 1:length(instr_forget_idx)
            forget_acts(:,i)=class_perf(:,instr_forget_idx(i));
    end
        
for_face = reshape(forget_acts(1,:),7,30)';
forget.acts.face = for_face;
for_scene =reshape(forget_acts(2,:),7,30)';
forget.acts.scene = for_scene;
for_object =reshape(forget_acts(3,:),7,30)';
forget.acts.object = for_object;
for_word =reshape(forget_acts(4,:),7,30)';
forget.acts.word = for_word;
for_rest =   reshape(forget_acts(5,:),7,30)';
forget.acts.rest = for_rest;
for_check = reshape(instr_remember_idx,7,30)';
forget.acts.check = for_check;

    for i = 1:length(instr_remember_neu_idx)
            remember_neu_acts(:,i)=class_perf(:,instr_remember_neu_idx(i));
    end

rem_neu_face = reshape(remember_neu_acts(1,:),7,15)';
rem_neu.acts.face = rem_neu_face;
rem_neu_scene =reshape(remember_neu_acts(2,:),7,15)';
rem_neu.acts.scene = rem_neu_scene;
rem_neu_object =reshape(remember_neu_acts(3,:),7,15)';
rem_neu.acts.object = rem_neu_object;
rem_neu_word =reshape(remember_neu_acts(4,:),7,15)';
rem_neu.acts.word = rem_neu_word;
rem_neu_rest =   reshape(remember_neu_acts(5,:),7,15)';
rem_neu.acts.rest = rem_neu_rest;
rem_neu_check = reshape(instr_remember_neu_idx,7,15)';
rem_neu.acts.check = rem_neu_check;

    for i = 1:length(instr_remember_neg_idx)
            remember_neg_acts(:,i)=class_perf(:,instr_remember_neg_idx(i));
    end

rem_neg_face = reshape(remember_neg_acts(1,:),7,15)';
rem_neg.acts.face = rem_neg_face;
rem_neg_scene =reshape(remember_neg_acts(2,:),7,15)';
rem_neg.acts.scene = rem_neg_scene;
rem_neg_object =reshape(remember_neg_acts(3,:),7,15)';
rem_neg.acts.object = rem_neg_object;
rem_neg_word =reshape(remember_neg_acts(4,:),7,15)';
rem_neg.acts.word = rem_neg_word;
rem_neg_rest =   reshape(remember_neg_acts(5,:),7,15)';
rem_neg.acts.rest = rem_neg_rest;
rem_neg_check = reshape(instr_remember_neg_idx,7,15)';
rem_neg.acts.check = rem_neg_check;

    for i = 1:length(instr_forget_neu_idx)
            forget_neu_acts(:,i)=class_perf(:,instr_forget_neu_idx(i));
    end

forget_neu_face = reshape(forget_neu_acts(1,:),7,15)';
forget_neu.acts.face = forget_neu_face;
forget_neu_scene =reshape(forget_neu_acts(2,:),7,15)';
forget_neu.acts.scene = forget_neu_scene;
forget_neu_object =reshape(forget_neu_acts(3,:),7,15)';
forget_neu.acts.object = forget_neu_object;
forget_neu_word =reshape(forget_neu_acts(4,:),7,15)';
forget_neu.acts.word = forget_neu_word;
forget_neu_rest =   reshape(forget_neu_acts(5,:),7,15)';
forget_neu.acts.rest = forget_neu_rest;
forget_neu_check = reshape(instr_forget_neu_idx,7,15)';
forget_neu.acts.check = forget_neu_check;

    for i = 1:length(instr_forget_neg_idx)
            forget_neg_acts(:,i)=class_perf(:,instr_forget_neg_idx(i));
    end

forget_neg_face = reshape(forget_neg_acts(1,:),7,15)';
forget_neg.acts.face = forget_neg_face;
forget_neg_scene =reshape(forget_neg_acts(2,:),7,15)';
forget_neg.acts.scene = forget_neg_scene;
forget_neg_object =reshape(forget_neg_acts(3,:),7,15)';
forget_neg.acts.object = forget_neg_object;
forget_neg_word =reshape(forget_neg_acts(4,:),7,15)';
forget_neg.acts.word = forget_neg_word;
forget_neg_rest =   reshape(forget_neg_acts(5,:),7,15)';
forget_neg.acts.rest = forget_neg_rest;
forget_neg_check = reshape(instr_forget_neg_idx,7,15)';
forget_neg.acts.check = forget_neg_check;


%%% averaging %%%
            
rem.mean.face = mean(rem_face);
rem.mean.scene = mean(rem_scene);
rem.mean.object = mean(rem_object);
rem.mean.word = mean(rem_word);
rem.mean.rest = mean(rem_rest);
   
forget.mean.face = mean(for_face);
forget.mean.scene = mean(for_scene);
forget.mean.object = mean(for_object);
forget.mean.word = mean(for_word);
forget.mean.rest = mean(for_rest);

rem_neu.mean.face = mean(rem_neu_face);
rem_neu.mean.scene = mean(rem_neu_scene);
rem_neu.mean.object = mean(rem_neu_object);
rem_neu.mean.word = mean(rem_neu_word);
rem_neu.mean.rest = mean(rem_neu_rest);

rem_neg.mean.face = mean(rem_neg_face);
rem_neg.mean.scene = mean(rem_neg_scene);
rem_neg.mean.object = mean(rem_neg_object);
rem_neg.mean.word = mean(rem_neg_word);
rem_neg.mean.rest = mean(rem_neg_rest);

forget_neu.mean.face = mean(forget_neu_face);
forget_neu.mean.scene = mean(forget_neu_scene);
forget_neu.mean.object = mean(forget_neu_object);
forget_neu.mean.word = mean(forget_neu_word);
forget_neu.mean.rest = mean(forget_neu_rest);

forget_neg.mean.face = mean(forget_neg_face);
forget_neg.mean.scene = mean(forget_neg_scene);
forget_neg.mean.object = mean(forget_neg_object);
forget_neg.mean.word = mean(forget_neg_word);
forget_neg.mean.rest = mean(forget_neg_rest);

%%standard error%%

rem.std.face = std(rem_face);
rem.std.scene = std(rem_scene);
rem.std.object = std(rem_object);
rem.std.word = std(rem_word);
rem.std.rest = std(rem_rest);
   
forget.std.face = std(for_face);
forget.std.scene = std(for_scene);
forget.std.object = std(for_object);
forget.std.word = std(for_word);
forget.std.rest = std(for_rest);

rem_neu.std.face = std(rem_neu_face);
rem_neu.std.scene = std(rem_neu_scene);
rem_neu.std.object = std(rem_neu_object);
rem_neu.std.word = std(rem_neu_word);
rem_neu.std.rest = std(rem_neu_rest);

rem_neg.std.face = std(rem_neg_face);
rem_neg.std.scene = std(rem_neg_scene);
rem_neg.std.object = std(rem_neg_object);
rem_neg.std.word = std(rem_neg_word);
rem_neg.std.rest = std(rem_neg_rest);

forget_neu.std.face = std(forget_neu_face);
forget_neu.std.scene = std(forget_neu_scene);
forget_neu.std.object = std(forget_neu_object);
forget_neu.std.word = std(forget_neu_word);
forget_neu.std.rest = std(forget_neu_rest);

forget_neg.std.face = std(forget_neg_face);
forget_neg.std.scene = std(forget_neg_scene);
forget_neg.std.object = std(forget_neg_object);
forget_neg.std.word = std(forget_neg_word);
forget_neg.std.rest = std(forget_neg_rest);

results.DF.remember = rem;
results.DF.forget = forget;
results.DF.neutral.remember = rem_neu;
results.DF.neutral.forget = forget_neu;
results.DF.negative.remember = rem_neg;
results.DF.negative.forget = forget_neg;


filename = sprintf('results_%s.mat',test_phase);
save(filename,'results');



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