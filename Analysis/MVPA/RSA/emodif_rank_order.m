function [rank_results] = emodif_rank_order(subjNum)

%rank order analysis - takes the output of intraphase RSA and rank orders
%(to output percentile ranking) the correlation to matching trials vs
%non-matching trials. 

% first part of this analysis is subject specific - but can do 'aggregate'
% in place of subjNum

args.base_dir = '/Users/tw24955/emodif_data';


%single subject
if subjNum == 'aggregate'

load('/Users/tw24955/EmoDif/Analysis/MVPA/subj_list.mat')

else
    
for i = 1:length(subj_list);
subjNum = (subj_list{i,1});
subj_folder = subj_list{i,2};
args.subjID = sprintf('emodif_%s',num2str(subjNum));
args.maskName = maskName;
args.output_dir = sprintf('%s/%s/%s/rank_order', args.base_dir, args.subjID, maskName);


args.output_dir = sprintf('%s/aggregate_results/RSA_results/%s', args.base_dir, maskName);
args.outfname = sprintf('%s/rsa_rankorder_aggregate_%s', args.output_dir, date);
args.localizer.trialnum = 24; % 30 for subjects 1-3, but already coded in ...
%within phase rsa script. 24 should apply to everyone else. reducing to 24 happens at concatenation level. comment out 1-3 
%if not using these first 4 subjects.
args.preview.trialnum = 60;
args.DFencode.trialnum = 60;

mkdir(args.output_dir);

args.subj_dir = sprintf('%s/%s', args.base_dir, args.subjID);

if type == 1 | type == 2 | type == 3
    args.data_dir = sprintf('%s/results/rsa_results/intraphase/%s/%d/%s',args.subj_dir,args.maskName, TRstart, subj_folder);
    
elseif type == 4
    
    args.data_dir = sprintf('%s/results/rsa_results/preview_dfencode/%s/%d/%s',args.subj_dir,args.maskName, TRstart, subj_folder);
end


results.bysubject.name = sprintf('%s/emodif_%d_TR%dto%d_rsa_results_results.mat',args.data_dir, subjNum, TRstart, ((TRstart+TRlength)-1));
results.bysubject.data(i) = load(results.bysubject.name);
results.bysubject.names{i} = results.bysubject.name;




end

    
end





