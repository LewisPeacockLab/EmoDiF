function [] = emodif_create_aggregate_fsor(maskName)

%emodif_create_aggregate_fsor('tempoccfusi_pHg_LOC_combined_epi_space')

%load subject folder list

load('/Users/tw24955/emodif_data/subj_folders_fsor.mat');
subj_list = subj_folders;

%subj_list

for i = 1:length(subj_list)
args.subjnum = (subj_list{i,1});
subj_DFencode = (subj_list{i,2});

args.data_dir = '/Users/tw24955/emodif_data';

args.maskName = maskName;
args.script_dir = pwd;

args.DFencode.trialnum = 60;

args.subjID = sprintf('emodif_%d',args.subjnum);

args.subj_dir = sprintf('%s/%s', args.data_dir, args.subjID);

results.bysubject.name = sprintf('%s/results/DFencode/%s/%s/results_DFencode_bar_compare.mat',args.subj_dir, maskName, subj_DFencode);
results.bysubject.data(i) = load(results.bysubject.name);
results.bysubject.names{i} = results.bysubject.name;
end

outfilename1 = 'emodif_aggregate_results_fsor_complete.mat';


save(outfilename1,'results');
end
