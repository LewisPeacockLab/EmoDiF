function [] = emodif_create_aggregate(maskName)

%emodif_create_aggregate('tempoccfusi_pHg_LOC_combined_epi_space')

%load subject folder list

load('/Users/tw24955/emodif_data/subj_folders_011419.mat');
subj_list = subj_folders;

%subj_list

for i = 1:length(subj_list);
args.subjID = (subj_list{i,1});
subj_DFencode = (subj_list{i,2});

args.data_dir = '/Users/tw24955/emodif_data';

args.maskName = maskName;
args.script_dir = pwd;

args.DFencode.trialnum = 60;

args.subj_dir = sprintf('%s/%s', args.data_dir, args.subjID);

results.bysubject.name = sprintf('%s/results/DFencode/%s/%s/results_DFencode_bar_compare.mat',args.subj_dir, maskName, subj_DFencode);
results.bysubject.data(i) = load(results.bysubject.name);
results.bysubject.names{i} = results.bysubject.name;

%create scene-word separation score using the postcue (TR 5, 6 and 7) by
%trial. 

results.bysubject.data(i).results.forget.all.bar.postcue.separation_bytrial = results.bysubject.data(i).results.forget.all.bar.postcue.scene.bytrial - results.bysubject.data(i).results.forget.all.bar.postcue.word.bytrial;
results.bysubject.data(i).results.remember.all.bar.postcue.separation_bytrial = results.bysubject.data(i).results.remember.all.bar.postcue.scene.bytrial - results.bysubject.data(i).results.remember.all.bar.postcue.word.bytrial;

results_separation.bysubject(i).remember = results.bysubject.data(i).results.remember.all.bar.postcue.separation_bytrial;
results_separation.bysubject(i).forget = results.bysubject.data(i).results.forget.all.bar.postcue.separation_bytrial;
end

outfilename1 = 'emodif_aggregate_results_complete.mat';
outfilename2 = 'emodif_aggregate_results_separation.mat';

save(outfilename1,'results');
save(outfilename2,'results_separation');

end
