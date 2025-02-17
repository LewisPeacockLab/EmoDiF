function emodif_logreg_boot(num_boot)
%to where results directory is
load('emodif_aggregate_results_separation_beh.mat')

for i = 1:num_boot

%create logreg - set up matrix of 720 R and 720 F trials 
R_data = [];
F_data = [];
R_behav = [];
F_behav = [];

for x = 1:24
    r = randi(24);
    r_data = results_separation.bysubject(r).remember(:);
    f_data = results_separation.bysubject(r).forget(:);
    r_behav = results_separation.bysubject(r).remember_beh_acc(:);
    f_behav = results_separation.bysubject(r).forget_beh_acc(:);
    R_data = vertcat(r_data, R_data);
    F_data = vertcat(f_data, F_data);
    R_behav = vertcat(r_behav, R_behav);
    F_behav = vertcat(f_behav, F_behav);

end
    R_behav_cat = categorical(R_behav);
    F_behav_cat = categorical(F_behav);
    results(i).R_data = R_data;
    results(i).R_behav_cat = R_behav_cat;
    results(i).F_data = F_data;
    results(i).F_behav_cat = F_behav_cat;
    
    
    
    [results(i).R_results,results(i).Rdev,results(i).Rstats] = mnrfit(R_data,R_behav_cat);
    [results(i).F_results,results(i).Fdev,results(i).Fstats] = mnrfit(F_data, F_behav_cat);
    
    results(i).R_results_beta = results(i).R_results(2,:);
    results(i).F_results_beta = results(i).F_results(2,:);
    
end

outfilename = sprintf('emodif_logreg_boot_%d.mat', num_boot);
save(outfilename, 'results')
end
    

