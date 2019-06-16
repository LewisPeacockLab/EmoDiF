function emodif_logreg_boot_stats(bootnum)

load(sprintf('emodif_logreg_boot_%d.mat',bootnum))

R_results_beta = vertcat(results.R_results_beta);
F_results_beta = vertcat(results.F_results_beta);
results_agg.R_results_beta = R_results_beta;
results_agg.F_results_beta = F_results_beta;

RminusF_results_beta = R_results_beta - F_results_beta;
results_agg.RminusF_results_beta = RminusF_results_beta;


R_results_beta_mean = mean(R_results_beta);
F_results_beta_mean = mean(F_results_beta);
RminusF_results_beta_mean = mean(RminusF_results_beta);
result_agg.R_results_beta_mean = R_results_beta_mean;
result_agg.F_results_beta_mean= F_results_beta_mean;

R_index = R_results_beta > 0;
F_index = F_results_beta < 0;
RminusF_index = RminusF_results_beta > 0;

R_count = count(R_index);
F_count = count(F_index);
RminusF_count = count(RminusF_index);

R_prob = R_count/bootnum;
F_prob = F_count/bootnum;
RminusF_prob = RminusF_count/bootnum;
results_agg.R_prob = R_prob;
results_agg.F_prob = F_prob;
results_agg.RminusF_prob = RminusF_prob;

outfilename = sprintf('emodif_logreg_boot_%d_stats.mat',bootnum);
save(outfilename, 'results_agg')
end