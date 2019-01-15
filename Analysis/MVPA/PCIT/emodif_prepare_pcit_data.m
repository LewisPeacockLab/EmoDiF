function [pcit_data ] = emodif_prepare_pcit_data(target,data_folders,mask)
% e.g.  emodif_prepare_pcit_data('words','subj_folders_011419', 'tempoccfusi_pHg_LOC_combined_epi_space')
% this generates PCIT data for face, scenes, objects, word and rest. This
% also generates for each type of stimuli, 2 conditions: remember and
% forget. 
% Data format for P-CIT
% nTrials rows x 6 cols
% 1 - Subject ID
% 2 - Trial number
% 3 - Category * for this experiement this will be stimuli type
% 4 - Predictor variable
% 5 - Dependent variable or BEHAVIOR
% 6 - Net effect cluster pre and post cue TR1,2,3 and TR 5,6,7. 

par.version =  '2019Jan14';
%target = words, scenes, w-s, s-w 

%%bring in subject list 
par.data_dir = '/Users/tw24955/emodif_data/';
par.PCITdata_dir = '/Users/tw24955/emodif_data/PCIT/data/';
par.data_folders= load(sprintf('%s/%s',par.data_dir,data_folders));
subjids = par.data_folders.subj_folders(:,1);
data_folders = par.data_folders.subj_folders(:,2);
par.mask = mask;


PCIT_rem_face = [];
PCIT_rem_scene = [];
PCIT_rem_object = [];
PCIT_rem_word = [];
PCIT_rem_rest = [];
PCIT_for_face = [];
PCIT_for_scene = [];
PCIT_for_object = [];
PCIT_for_word = [];
PCIT_for_rest = [];

for i = 1:length(subjids)
    subj_id = subjids{i};
    data_folder = data_folders{i};

par.subj_dir = sprintf('%s/%s',par.data_dir, subj_id);
par.results_folder = sprintf('%s/results/DFencode/%s/%s',par.subj_dir,par.mask, data_folder);
load(sprintf('%s/results_DFencode.mat',par.results_folder));

%convert behavior into scale from 0 to 1

for b = 1:length(results.remember.all.beh)
    if results.remember.all.beh(b) == 4
        results.remember.all.behconvert(b) = 0 ;
    elseif results.remember.all.beh(b) == 3
        results.remember.all.behconvert(b) = .33;
    elseif results.remember.all.beh(b) == 2
        results.remember.all.behconvert(b) = .66;
    elseif results.remember.all.beh(b) == 1
        results.remember.all.behconvert(b) = 1;
    end
end

results.remember.all.behconvert = results.remember.all.behconvert';

for b = 1:length(results.forget.all.beh)
    if results.forget.all.beh(b) == 4
        results.forget.all.behconvert(b) = 0 ;
    elseif results.forget.all.beh(b) == 3
        results.forget.all.behconvert(b) = .33;
    elseif results.forget.all.beh(b) == 2
        results.forget.all.behconvert(b) = .66;
    elseif results.forget.all.beh(b) == 1
        results.forget.all.behconvert(b) = 1;
    end
end

results.forget.all.behconvert = results.forget.all.behconvert';

%build PCIT data structure for each condition REMEBER and FORGET with
%conditions 1:Face 2:Scene 3:Object 4:word 5:Rest
%each subject will have 2 PCIT data structures for aggregation

face_rem = results.remember.all.acts.face;
face_for = results.forget.all.acts.face;

face_rem_1 =  face_rem(:,1:3);
face_rem_pre = mean(face_rem_1,2);

face_rem_2 = face_rem(:,5:7);
face_rem_post = mean(face_rem_2,2);

face_for_1 =  face_for(:,1:3);
face_for_pre = mean(face_for_1,2);

face_for_2 = face_for(:,5:7);
face_for_post = mean(face_for_2,2);


scene_rem = results.remember.all.acts.scene;
scene_for = results.forget.all.acts.scene;

scene_rem_1 =  scene_rem(:,1:3);
scene_rem_pre = mean(scene_rem_1,2);

scene_rem_2 = scene_rem(:,5:7);
scene_rem_post = mean(scene_rem_2,2);

scene_for_1 =  scene_for(:,1:3);
scene_for_pre = mean(scene_for_1,2);

scene_for_2 = scene_for(:,5:7);
scene_for_post = mean(scene_for_2,2);


object_rem = results.remember.all.acts.object;
object_for = results.forget.all.acts.object;

object_rem_1 =  object_rem(:,1:3);
object_rem_pre = mean(object_rem_1,2);

object_rem_2 = object_rem(:,5:7);
object_rem_post = mean(object_rem_2,2);

object_for_1 =  object_for(:,1:3);
object_for_pre = mean(object_for_1,2);

object_for_2 = object_for(:,5:7);
object_for_post = mean(object_for_2,2);

word_rem = results.remember.all.acts.word;
word_for = results.forget.all.acts.word;

word_rem_1 =  word_rem(:,1:3);
word_rem_pre = mean(word_rem_1,2);

word_rem_2 = word_rem(:,5:7);
word_rem_post = mean(word_rem_2,2);

word_for_1 =  word_for(:,1:3);
word_for_pre = mean(word_for_1,2);

word_for_2 = word_for(:,5:7);
word_for_post = mean(word_for_2,2);

rest_rem = results.remember.all.acts.rest;
rest_for = results.forget.all.acts.rest;

rest_rem_1 =  rest_rem(:,1:3);
rest_rem_pre = mean(rest_rem_1,2);

rest_rem_2 = rest_rem(:,5:7);
rest_rem_post = mean(rest_rem_2,2);

rest_for_1 =  rest_for(:,1:3);
rest_for_pre = mean(rest_for_1,2);

rest_for_2 = rest_for(:,5:7);
rest_for_post = mean(rest_for_2,2);

% now build each trial pre and post effect will end up with 10 matrices - 
% R - 1. face, 2. scene, 3. object, 4. word, 5. rest
% F - 6. face, 7. scene, 8. object, 9. word, 10. rest

PCIT_rem_face_sub = [];
PCIT_rem_scene_sub = [];
PCIT_rem_object_sub = [];
PCIT_rem_word_sub = [];
PCIT_rem_rest_sub = [];
PCIT_for_face_sub = [];
PCIT_for_scene_sub = [];
PCIT_for_object_sub = [];
PCIT_for_word_sub = [];
PCIT_for_rest_sub = [];

for t = 1:30
    PCIT_rem_face_s(1,1) = i; % first column is subject
    PCIT_rem_face_s(1,2) = t;
    PCIT_rem_face_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_face_s(1,4) = face_rem_pre(t); % net pre
    PCIT_rem_face_s(1,5) = results.remember.all.behconvert(t);
    PCIT_rem_face_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_face_s(2,1) = i;
    PCIT_rem_face_s(2,2) = t;
    PCIT_rem_face_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_face_s(2,4) = face_rem_post(t); % net pre
    PCIT_rem_face_s(2,5) = results.remember.all.behconvert(t);
    PCIT_rem_face_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_face_sub = vertcat(PCIT_rem_face_sub, PCIT_rem_face_s);
end

for t = 1:30
    PCIT_rem_scene_s(1,1) = i; % first column is subject
    PCIT_rem_scene_s(1,2) = t;
    PCIT_rem_scene_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_scene_s(1,4) = scene_rem_pre(t); % net pre
    PCIT_rem_scene_s(1,5) = results.remember.all.behconvert(t);
    PCIT_rem_scene_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_scene_s(2,1) = i;
    PCIT_rem_scene_s(2,2) = t;
    PCIT_rem_scene_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_scene_s(2,4) = scene_rem_post(t); % net pre
    PCIT_rem_scene_s(2,5) = results.remember.all.behconvert(t);
    PCIT_rem_scene_s(2,6) = t+(30*(i-1));%net effects, same as trial number
    
    PCIT_rem_scene_sub = vertcat(PCIT_rem_scene_sub, PCIT_rem_scene_s);
end

for t = 1:30
    PCIT_rem_object_s(1,1) = i; % first column is subject
    PCIT_rem_object_s(1,2) = t;
    PCIT_rem_object_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_object_s(1,4) = object_rem_pre(t); % net pre
    PCIT_rem_object_s(1,5) = results.remember.all.behconvert(t);
    PCIT_rem_object_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_object_s(2,1) = i;
    PCIT_rem_object_s(2,2) = t;
    PCIT_rem_object_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_object_s(2,4) = object_rem_post(t); % net pre
    PCIT_rem_object_s(2,5) = results.remember.all.behconvert(t);
    PCIT_rem_object_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_object_sub = vertcat(PCIT_rem_object_sub, PCIT_rem_object_s);
end

for t = 1:30
    PCIT_rem_word_s(1,1) = i; % first column is subject
    PCIT_rem_word_s(1,2) = t;
    PCIT_rem_word_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_word_s(1,4) = word_rem_pre(t); % net pre
    PCIT_rem_word_s(1,5) = results.remember.all.behconvert(t);
    PCIT_rem_word_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_word_s(2,1) = i;
    PCIT_rem_word_s(2,2) = t;
    PCIT_rem_word_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_word_s(2,4) = word_rem_post(t); % net pre
    PCIT_rem_word_s(2,5) = results.remember.all.behconvert(t);
    PCIT_rem_word_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_word_sub = vertcat(PCIT_rem_word_sub, PCIT_rem_word_s);
end

for t = 1:30
    PCIT_rem_rest_s(1,1) = i; % first column is subject
    PCIT_rem_rest_s(1,2) = t;
    PCIT_rem_rest_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_rest_s(1,4) = rest_rem_pre(t); % net pre
    PCIT_rem_rest_s(1,5) = results.remember.all.behconvert(t);
    PCIT_rem_rest_s(1,6) = t+(30*(i-1));%net effects, same as trial number
    
    PCIT_rem_rest_s(2,1) = i;
    PCIT_rem_rest_s(2,2) = t;
    PCIT_rem_rest_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_rem_rest_s(2,4) = rest_rem_post(t); % net pre
    PCIT_rem_rest_s(2,5) = results.remember.all.behconvert(t);
    PCIT_rem_rest_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_rem_rest_sub = vertcat(PCIT_rem_rest_sub, PCIT_rem_rest_s);
end

for t = 1:30
    PCIT_for_face_s(1,1) = i; % first column is subject
    PCIT_for_face_s(1,2) = t;
    PCIT_for_face_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_for_face_s(1,4) = face_for_pre(t); % net pre
    PCIT_for_face_s(1,5) = results.forget.all.behconvert(t);
    PCIT_for_face_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_face_s(2,1) = i;
    PCIT_for_face_s(2,2) = t;
    PCIT_for_face_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_for_face_s(2,4) = face_for_post(t); % net pre
    PCIT_for_face_s(2,5) = results.forget.all.behconvert(t);
    PCIT_for_face_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_face_sub = vertcat(PCIT_for_face_sub, PCIT_for_face_s);
end

for t = 1:30
    PCIT_for_scene_s(1,1) = i; % first column is subject
    PCIT_for_scene_s(1,2) = t;
    PCIT_for_scene_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_for_scene_s(1,4) = scene_for_pre(t); % net pre
    PCIT_for_scene_s(1,5) = results.forget.all.behconvert(t);
    PCIT_for_scene_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_scene_s(2,1) = i;
    PCIT_for_scene_s(2,2) = t;
    PCIT_for_scene_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_for_scene_s(2,4) = scene_for_post(t); % net pre
    PCIT_for_scene_s(2,5) = results.forget.all.behconvert(t);
    PCIT_for_scene_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_scene_sub = vertcat(PCIT_for_scene_sub, PCIT_for_scene_s);
end

for t = 1:30
    PCIT_for_object_s(1,1) = i; % first column is subject
    PCIT_for_object_s(1,2) = t;
    PCIT_for_object_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_for_object_s(1,4) = object_for_pre(t); % net pre
    PCIT_for_object_s(1,5) = results.forget.all.behconvert(t);
    PCIT_for_object_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_object_s(2,1) = i;
    PCIT_for_object_s(2,2) = t;
    PCIT_for_object_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_for_object_s(2,4) = object_for_post(t); % net pre
    PCIT_for_object_s(2,5) = results.forget.all.behconvert(t);
    PCIT_for_object_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_object_sub = vertcat(PCIT_for_object_sub, PCIT_for_object_s);
end

for t = 1:30
    PCIT_for_word_s(1,1) = i; % first column is subject
    PCIT_for_word_s(1,2) = t;
    PCIT_for_word_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_for_word_s(1,4) = word_for_pre(t); % net pre
    PCIT_for_word_s(1,5) = results.forget.all.behconvert(t);
    PCIT_for_word_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_word_s(2,1) = i;
    PCIT_for_word_s(2,2) = t;
    PCIT_for_word_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_for_word_s(2,4) = word_for_post(t); % net pre
    PCIT_for_word_s(2,5) = results.forget.all.behconvert(t);
    PCIT_for_word_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_word_sub = vertcat(PCIT_for_word_sub, PCIT_for_word_s);
end

for t = 1:30
    PCIT_for_rest_s(1,1) = i; % first column is subject
    PCIT_for_rest_s(1,2) = t;
    PCIT_for_rest_s(1,3) = 1; % there is only one category for these PCIT data
    PCIT_for_rest_s(1,4) = rest_for_pre(t); % net pre
    PCIT_for_rest_s(1,5) = results.forget.all.behconvert(t);
    PCIT_for_rest_s(1,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_rest_s(2,1) = i;
    PCIT_for_rest_s(2,2) = t;
    PCIT_for_rest_s(2,3) = 1; % there is only one category for these PCIT data
    PCIT_for_rest_s(2,4) = rest_for_post(t); % net pre
    PCIT_for_rest_s(2,5) = results.forget.all.behconvert(t);
    PCIT_for_rest_s(2,6) = t+(30*(i-1)); %net effects, same as trial number
    
    PCIT_for_rest_sub = vertcat(PCIT_for_rest_sub, PCIT_for_rest_s);
end

PCIT_rem_face = vertcat(PCIT_rem_face, PCIT_rem_face_sub);
PCIT_rem_scene = vertcat(PCIT_rem_scene, PCIT_rem_scene_sub);
PCIT_rem_object = vertcat(PCIT_rem_object, PCIT_rem_object_sub);
PCIT_rem_word = vertcat(PCIT_rem_word, PCIT_rem_word_sub);
PCIT_rem_rest = vertcat(PCIT_rem_rest, PCIT_rem_rest_sub);
PCIT_for_face = vertcat(PCIT_for_face, PCIT_for_face_sub);
PCIT_for_scene = vertcat(PCIT_for_scene, PCIT_for_scene_sub);
PCIT_for_object = vertcat(PCIT_for_object, PCIT_for_object_sub);
PCIT_for_word = vertcat(PCIT_for_word, PCIT_for_word_sub);
PCIT_for_rest = vertcat(PCIT_for_rest, PCIT_for_rest_sub);

end

cd(par.PCITdata_dir)

outfilename_rem_face = 'emodif_PCIT_face_remember_pre_3TRs_post_3TRs.mat';
data = PCIT_rem_face;
save(outfilename_rem_face, 'data');

outfilename_rem_scene = 'emodif_PCIT_scene_remember_pre_3TRs_post_3TRs.mat';
data = PCIT_rem_scene;
save(outfilename_rem_scene, 'data');

outfilename_rem_object = 'emodif_PCIT_object_remember_pre_3TRs_post_3TRs.mat';
data = PCIT_rem_object;
save(outfilename_rem_object, 'data');

outfilename_rem_word  = 'emodif_PCIT_word_remember_pre_3TRs_post_3TRs.mat';
data = PCIT_rem_word;
save(outfilename_rem_word, 'data');

outfilename_rem_rest  = 'emodif_PCIT_rest_remember_pre_3TRs_post_3TRs.mat';
data = PCIT_rem_rest;
save(outfilename_rem_rest, 'data');

outfilename_for_face = 'emodif_PCIT_face_forget_pre_3TRs_post_3TRs.mat';
data = PCIT_for_face;
save(outfilename_for_face, 'data');

outfilename_for_scene = 'emodif_PCIT_scene_forget_pre_3TRs_post_3TRs.mat';
data = PCIT_for_scene;
save(outfilename_for_scene, 'data');

outfilename_for_object = 'emodif_PCIT_object_forget_pre_3TRs_post_3TRs.mat';
data = PCIT_for_object;
save(outfilename_for_object, 'data');

outfilename_for_word  = 'emodif_PCIT_word_forget_pre_3TRs_post_3TRs.mat';
data = PCIT_for_word;
save(outfilename_for_word, 'data');

outfilename_for_rest  = 'emodif_PCIT_rest_forget_pre_3TRs_post_3TRs.mat';
data = PCIT_for_rest;
save(outfilename_for_rest, 'data');


end









