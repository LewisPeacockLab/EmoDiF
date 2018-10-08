function [rank_results] = emodif_rank_order(maskName, aggregate_file_date, preview_shift, TRstart, TRlength)

% [rank_results] = = emodif_rank_order('scene-all_005', '02-Oct-2018', 2, 2, 3)
%rank order analysis - takes the output of intraphase RSA and rank orders
%(to output percentile ranking) the correlation to matching trials vs
%non-matching trials. 

% first part of this analysis is subject specific - but can do 'aggregate'
% in place of subjNum

args.base_dir = '/Users/tw24955/emodif_data';


%single subject
% if subjNum == 'aggregate'
% 
load(sprintf('/Users/tw24955/EmoDif/Analysis/MVPA/subj_list_%s.mat', maskName));
% 
% else

args.output_dir = sprintf('%s/aggregate_results/RSA_rank_results/%s', args.base_dir, maskName);
args.outfname = sprintf('%s/rsa_rankorder_aggregate_%s', args.output_dir, date);
args.localizer.trialnum = 24; % 30 for subjects 1-3, but already coded in ...
%within phase rsa script. 24 should apply to everyone else. reducing to 24 happens at concatenation level. comment out 1-3 
%if not using these first 4 subjects.
args.preview.trialnum = 60;
args.DFencode.trialnum = 60;

args.data_dir = sprintf('%s/aggregate_results/RSA_preview_dfencode_results/%s',args.base_dir,maskName);

mkdir(args.output_dir);
args.maskName = maskName;
args.aggregate_file_date = aggregate_file_date;
args.data_file = sprintf('%s/rsa_preview_DFencode_aggregate_%s',args.data_dir, args.aggregate_file_date);

load(args.data_file);

for i = 1:length(results.bysubject.names);
corr_matrix_match_fullz = results.bysubject.data(i).rsa.results.smatrix.corr_matrix_match_fullz;
corr_matrix_match_Fz = results.bysubject.data(i).rsa.results.smatrix.corr_matrix_match_Fz;
corr_matrix_match_Rz = results.bysubject.data(i).rsa.results.smatrix.corr_matrix_match_Rz;

for x = 1:length(corr_matrix_match_fullz(1,:));
    order = (1:length(corr_matrix_match_fullz(1,:)))';
    rank = (1:length(corr_matrix_match_fullz(1,:)))';
    to_sort = horzcat(corr_matrix_match_fullz(:,x),order);
    sorted_corr_matrix_match_fullz_x = sortrows(to_sort);
    ranked_corr_matrix_match_fullz_x = horzcat(rank, sorted_corr_matrix_match_fullz_x); % rank, correlation, order
    match = find(ranked_corr_matrix_match_fullz_x (:,3) == x);
    matched_rank = ranked_corr_matrix_match_fullz_x(match,1)/length(corr_matrix_match_fullz(1,:));
    results.rankorder.rank_corr_matrix_match_fullz(i,x) = matched_rank;
end
    
for x = 1:length(corr_matrix_match_Fz(1,:));
    order = (1:length(corr_matrix_match_Fz(1,:)))';
    rank = (1:length(corr_matrix_match_Fz(1,:)))';
    to_sort = horzcat(corr_matrix_match_Fz(:,x),order);
    sorted_corr_matrix_match_Fz_x = sortrows(to_sort);
    ranked_corr_matrix_match_Fz_x = horzcat(rank, sorted_corr_matrix_match_Fz_x); % rank, correlation, order
    match = find(ranked_corr_matrix_match_Fz_x (:,3) == x);
    matched_rank = ranked_corr_matrix_match_Fz_x(match,1)/length(corr_matrix_match_Fz(1,:));
    results.rankorder.rank_corr_matrix_match_Fz(i,x) = matched_rank;
end

for x = 1:length(corr_matrix_match_Rz(1,:));
    order = (1:length(corr_matrix_match_Rz(1,:)))';
    rank = (1:length(corr_matrix_match_Rz(1,:)))';
    to_sort = horzcat(corr_matrix_match_Rz(:,x),order);
    sorted_corr_matrix_match_Rz_x = sortrows(to_sort);
    ranked_corr_matrix_match_Rz_x = horzcat(rank, sorted_corr_matrix_match_Rz_x); % rank, correlation, order
    match = find(ranked_corr_matrix_match_Rz_x (:,3) == x);
    matched_rank = ranked_corr_matrix_match_Rz_x(match,1)/length(corr_matrix_match_Rz(1,:));
    results.rankorder.rank_corr_matrix_match_Rz(i,x) = matched_rank;
end

end

%aggregate all subjects 


results.rankorder.rank_corr_matrix_match_fullz_mean_bysub = mean(results.rankorder.rank_corr_matrix_match_fullz,2);
results.rankorder.rank_corr_matrix_match_Rz_mean_bysub = mean(results.rankorder.rank_corr_matrix_match_Rz,2);
results.rankorder.rank_corr_matrix_match_Fz_mean_bysub = mean(results.rankorder.rank_corr_matrix_match_Fz,2);
results.rankorder.rank_corr_matrix_match_fullz_mean = mean(results.rankorder.rank_corr_matrix_match_fullz_mean_bysub);
results.rankorder.rank_corr_matrix_match_Rz_mean = mean(results.rankorder.rank_corr_matrix_match_Rz_mean_bysub);
results.rankorder.rank_corr_matrix_match_Fz_mean = mean(results.rankorder.rank_corr_matrix_match_Fz_mean_bysub);

filename = args.outfname;
save(filename, 'results');
end

    






