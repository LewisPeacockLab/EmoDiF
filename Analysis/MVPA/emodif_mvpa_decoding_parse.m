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
   args.subj_dir = sprintf('/Users/tw24955/emodif_data/%s', args.subjID);
% for tigger
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
    
%     for k = 1:length(args.DF.trialn)
%         remember_reg.face(1,(k-1):() 

    %this is where we should put the rest of the regressors as a sanity
    %cut regressors for 1 ST 
    %check
    for i = 1:length(instr_remember_idx)
            remember_acts(:,i)=class_perf(:,instr_remember_idx(i));
    end
        
rem_face = reshape(remember_acts(1,:),7,30)';
rem_scene =reshape(remember_acts(2,:),7,30)';
rem_object =reshape(remember_acts(3,:),7,30)';
rem_word =reshape(remember_acts(4,:),7,30)';
rem_rest =   reshape(remember_acts(5,:),7,30)';


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
results.DF.remember = rem;
results.DF.forget = forget;

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