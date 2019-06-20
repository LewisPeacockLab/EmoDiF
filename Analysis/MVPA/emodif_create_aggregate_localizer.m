function [] = emodif_create_aggregate_localizer(maskName)
%emodif_create_aggregate_localizer('tempoccfusi_pHg_LOC_combined_epi_space')

load('/Users/tw24955/emodif_data/subj_folders_011419.mat');
subj_list = subj_folders;

for i = 1:length(subj_list);
args.subjID = (subj_list{i,1});
args.localizer = '25-Sep-2018';

args.data_dir = '/Users/tw24955/emodif_data';

args.maskName = maskName;
args.script_dir = pwd;


args.subj_dir = sprintf('%s/%s', args.data_dir, args.subjID);

results.bysubject.name = sprintf('%s/results/localizer/%s/%s/localizer_parsed_corr.mat',args.subj_dir, maskName, args.localizer);
results.bysubject.data(i) = load(results.bysubject.name);
results.bysubject.names{i} = results.bysubject.name;

end

outfilename1 = 'emodif_aggregate_localizer_corr.mat';
save(outfilename1,'results');